library(MASS)
library(Rcpp)
library(geepack)
source("/Users/lanluo/Google Drive/RenewEE/R code/QIF_cpp.r")
source("/Users/lanluo/Google Drive/RenewEE/R code/Project 3/data_cs.r")

sourceCpp("/Users/lanluo/Google Drive/RenewEE/R code/Project 3/increQIF_cs.cpp")
source("/Users/lanluo/Google Drive/RenewEE/R code/Project 3/RenewQIF_cs.R")

#sourceCpp("/Users/lanluo/Google Drive/RenewEE/R code/increQIF.cpp")
sourceCpp("/Users/lanluo/Google Drive/RenewEE/R code/increQIF_sub.cpp")



##########################################
## Evaluation functions
##########################################
eval.func <- function(beta, betahat, sd){
    
    bias<-abs(betahat - beta)[2]
    sd <- sd[2]
    lower <- betahat[2] - 1.96 * sd
    upper <- betahat[2] + 1.96 * sd
    if(beta[2]>= lower && beta[2] <= upper){
    	covind <- 1
    }else{
    	covind <- 0
    }

    c(bias=bias, sd=sd, covind=covind)
}


setwd("/Users/lanluo/Google Drive/RenewEE/R code/Project 3")

#number of simulations
nsim<-500

# MB visits for each individual, a total of n individuals
m <- 5
B <- 50
MB <- m * B

p <- 3
n <- 500 

hete_rho<- FALSE

beta0<-c(0.2,-0.2,0.2)

type<- family <-"gaussian"
corst_x<-"cs"
rho_x <- 0.5

## True correlation matrix 
corst_y <- "cs"
rho_y <- 0.2

intercept <- TRUE
## Working matrix
corstr <- "exchangeable"
outlier <- FALSE
B_c <- 0
n_c <- 0

outputfilename <- paste(corst_y, rho_y, corstr, "m", m, "n", n, "B", B, "p", p, sep="")

for(s in c(1:nsim)){
    print(s)

    tempdatadir<-paste("Temp_",outputfilename, "_", s, sep="")
    datagenerator_cs(n, m, p, B, tempdatadir, type, beta0, intercept=TRUE, corst_x=corst_x, rho_x=rho_x, 
 corst_y=corst_y, rho_y=rho_y, seed=s)
  

    time.all.read<-system.time( load(paste(tempdatadir, "/full", ".RData", sep="")) )[3]

    # 1. GEE on full data
#    time.all.gee<-system.time(result.A.gee<-geeglm(y_com ~ X_com, id=id_com, family=family, corstr="ar1"))[3]+ time.all.read
#    print("GEE on full data: done!")

    # 2. QIF on full data
    time.all.qif <- system.time(result.A.qif <- QIF_cpp(X=X_com, y=y_com, id=id_com, 
      family=family, intercept=intercept, corstr=corstr))[3]+ time.all.read
    print("QIF on full data: done!")
    rm(y_com)
    rm(X_com)
    rm(id_com)
    
    # 3. QIF on sequential data batches
#    result.B.renewqif_ind <- renewqif_ind(B=B,tempdatadir=tempdatadir, family=family, intercept=TRUE, corstr=corstr)
#    time.B.renewqif_ind <- result.B.renewqif_ind$time
#    run.B.renewqif_ind <- result.B.renewqif_ind$run
#    print("QIF ind on sequential data: done!")

    result.B.renewqif_cs <- renewqif_cs(B=B,tempdatadir=tempdatadir, family=family, intercept=TRUE, corstr=corstr)
    time.B.renewqif_cs <- result.B.renewqif_cs$time
    run.B.renewqif_cs <- result.B.renewqif_cs$run
    print("QIF cs on sequential data: done!")


    unlink(tempdatadir, recursive = TRUE)

    out<-c(
    #    eval.func(beta0,coef(summary(result.A.gee))[,1],coef(summary(result.A.gee))[,2]), time.all.gee, time.all.gee-time.all.read,
        
        eval.func(beta0, result.A.qif$beta[,1], result.A.qif$beta[,2]), time.all.qif, time.all.qif-time.all.read,
    #    eval.func(beta0,result.B.renewqif_ind$beta[,1], result.B.renewqif_ind$beta[,2]),time.B.renewqif_ind, run.B.renewqif_ind,

        eval.func(beta0,result.B.renewqif_cs$beta[,1], result.B.renewqif_cs$beta[,2]),time.B.renewqif_cs, run.B.renewqif_cs

    )

   out_ese<-c(#coef(summary(result.A.gee))[,1],
   	result.A.qif$beta[,"Estimate"],
  #  result.B.renewqif_ind$beta[,"Estimate"], 
    result.B.renewqif_cs$beta[,"Estimate"]
    )
    
   # out_Q <- result.B.renewqif$Q

    write.table(as.matrix(t(out)), file=paste(outputfilename, ".csv", sep=""), sep=",",
        col.names=FALSE, row.names=s, append=TRUE)
    write.table(as.matrix(t(out_ese)), file=paste(outputfilename, "ese.csv", sep=""), sep=",",
      col.names=FALSE, row.names=s, append=TRUE)
   # write.table(as.matrix(t(out_Q)), file=paste(outputfilename,"Q.csv", sep=""), sep=",",
   #     col.names=FALSE, row.names=s, append=TRUE)
}



out=read.csv(paste(outputfilename, ".csv", sep=""), header=FALSE, col.names=c("num","gee_bias","gee_sd","gee_cov","gee_time","gee_run",
            "QIF_bias","QIF_sd","QIF_cov","QIF_time","QIF_run",
            "renewQIF_bias","renewQIF_sd","renewQIF_cov","renewQIF_time","renewQIF_run"))
colMeans(out)


estimate<-read.table(paste(outputfilename, "ese.csv", sep=""), sep=",")
col_sd=apply(estimate,2,sd)
gee_ese<-mean(col_sd[2:(2 + p-1)])
qif_ese<-mean(col_sd[(2 + p): (2 + p + p-1)])
renewqif_ese<-mean(col_sd[(2+p+p):(2+p+p+p-1)])
renewqifar1_ese<-mean(col_sd[(2+3*p):(2+3*p+p-1)])
c(gee_ese,qif_ese, renewqif_ese, renewqifar1_ese)


Q_homo <-read.table(paste("HomoAR-1m10N2000B40n50p2_gaussian_AR-10.5", "Q.csv", sep=""), sep=",")
# calculate column mean
col_Q_homo <- apply(Q_homo,2,mean)
pvalue_homo <- pchisq(col_Q_homo[-1],df=p,lower.tail=FALSE)

Q_hete <-read.table(paste(outputfilename, "Q.csv", sep=""), sep=",")
# calculate column mean
col_Q_hete <- apply(Q_hete,2,mean)
pvalue_hete <- pchisq(col_Q_hete[-1],df=p,lower.tail=FALSE)

par(mfrow=c(1,2))
par(mar=c(4.5,4.5,2.5,0.5))

plot(x=seq(1,B,1),y=result.B.renewgee$trace,type="b",xlab="Data batches", ylab="RenewGEE",lwd=1)
abline(v=3,col="red")
abline(h=0.5,col="blue",lwd=2)
plot(x=seq(1,B,1),y=result.B.renewqif$trace,type="b",xlab="Data batches", ylab="RenewQIF",lwd=1, ylim=c(0.5,1))
abline(v=3,col="red")
abline(h=0.5,col="blue",lwd=2)



