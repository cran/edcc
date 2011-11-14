.ecoxbar <- function(h=seq(0.1,1,by=.01),L=seq(2,4.5,by=.01),n=1:15,lambda=.05,delta=2,P0=NULL,P1=NULL,C0=NULL,C1=NULL,Cr=25,Cf=50,T0=0.0167,Tc=1,Tf=0,Tr=0,a=1,b=.1,d1=1,d2=1,call.print=TRUE,...){
  f <- function(h,L,n){
    alpha <- 2*pnorm(-L)
    beta <- pnorm(L-delta*sqrt(n))-pnorm(-L-delta*sqrt(n))
    ARL1 <- 1/alpha
    ARL2 <- 1/(1-beta)
    tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
    s <- exp(-lambda*h)/(1-exp(-lambda*h))
    if(!is.null(P0)&!is.null(P1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECP <- P0/lambda + P1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) - s*Cf/ARL1 - Cr - (a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      cost <- P0 - ECP/ECT
    }else
    if(!is.null(C0)&!is.null(C1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECC <- C0/lambda + C1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) + s*Cf/ARL1+Cr+(a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      cost <- ECC/ECT
    }else
    stop("You should at least give a pair of value to P0,P1 or C0,C1")
    return(cost)
  } 
  
  mat <- outer(h,L,FUN=f,n=n)
  aa <- which(mat==min(mat),arr.ind=TRUE)
  optimum <- c(n,L[aa[1,2]],h[aa[1,1]],min(mat))
  names(optimum) <- c("n","Optimum L","Optimum h","Cost")

    par(mar=c(7.1,4.1,2.1,2.1))
    contour(h,L,mat,xlab="h",ylab="L",...)
    points(optimum[3],optimum[2],pch=3)
    mtext(sprintf('n=%s   Opt L=%s   Opt h=%s   Cost=%s',optimum[1], optimum[2], optimum[3],round(optimum[4],digits=4)),side=1,line=4.5)
  if(call.print){
    ca <- match.call()
    ca[[1]] <- NULL
    mtext(paste("ecoXbar(",paste(names(ca),"=",unlist(ca),sep="",collapse=", "),")"),side=3,line=0)
  }

}



.ecocusum <- function( h=seq(.1,2,by=.1), H=seq(.1,1,by=.01),n=3:7,delta = 2,lambda = .01, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 20, Cf = 10,T0 = 0, Tc = .1,Tf = .1,Tr = 0.2, a = .5, b = .1,d1=0,d2=0,sided = "one",call.print=TRUE,...){
    f <- function(h,H,n){
    delta.std <- sqrt(n)*delta          #standardization for delta
    k <- delta.std/2
    ARL1 <- as.numeric(xcusum.arl(k,H,0,sided=sided))
    ARL2 <- as.numeric(xcusum.arl(k,H,delta.std,sided=sided))
    tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
    s <- exp(-lambda*h)/(1-exp(-lambda*h))
    if(!is.null(P0)&!is.null(P1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECP <- P0/lambda + P1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) - s*Cf/ARL1 - Cr - (a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      cost <- P0 - ECP/ECT
    }else
    if(!is.null(C0)&!is.null(C1)){
      ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
      ECC <- C0/lambda + C1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) + s*Cf/ARL1+Cr+(a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
      cost <- ECC/ECT
    }else
    stop("You should at least give a pair of value to P0,P1 or C0,C1")
    return(cost)
  }
      

    mat <- matrix(NA,length(h),length(H))
    for(i in 1:length(h))
      for(j in 1:length(H))
        mat[i,j] <- f(h[i],H[j],n)
    aa <- which(mat==min(mat),arr.ind=TRUE)
    optimum <- c(n,H[aa[1,][2]],h[aa[1,][1]],min(mat))
    names(optimum) <- c("n","Optimum H","Optimum h","Cost")

    par(mar=c(7.1,4.1,2.1,2.1))
    contour(h,H,mat,xlab="h",ylab="H",...)
    points(optimum[3],optimum[2],pch=3)
    mtext(sprintf('n=%s   Opt H=%s   Opt h=%s   Cost=%s',optimum[1], optimum[2], optimum[3],round(optimum[4],digits=4)),side=1,line=4.5)
    if(call.print){
      ca <- match.call()
      ca[[1]] <- NULL
      mtext(paste("ecoEwma(",paste(names(ca),"=",unlist(ca),sep="",collapse=", "),")"),side=3,line=0)
    }
  }

#.ecocusum(P0=150,P1=50,Cr=30,n=5,nlevels=50)

 .ecoewma <- function( h=seq(.1,2,by=.1), w=seq(0.01,1,by=.01),k=seq(2,4,by=0.1),n=3:8,delta = 2,lambda = .01, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 20, Cf = 10,T0 = 0,Tc = .1, Tf = .1, Tr = 0.2, a = .5, b = .1,d1=1,d2=1,sided="two",contour.plot=FALSE,call.print=TRUE,...){
   f <- function(h,w,k,n){
     delta.std <- sqrt(n)*delta #standardization fordelta
     ARL1 <- as.numeric(xewma.arl(w,k,0,sided=sided))
     ARL2 <- as.numeric(xewma.arl(w,k,delta.std,sided=sided))
     tau <- (1-(1+lambda*h)*exp(-lambda*h))/(lambda*(1-exp(-lambda*h)))
     s <- exp(-lambda*h)/(1-exp(-lambda*h))
     if(!is.null(P0)&!is.null(P1)){
       ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
       ECP <- P0/lambda + P1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) - s*Cf/ARL1 - Cr - (a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
       cost <- P0 - ECP/ECT
     }else
     if(!is.null(C0)&!is.null(C1)){
       ECT <- 1/lambda+(1-d1)*s*Tf/ARL1 - tau + n*T0 + h*ARL2 + Tc + Tr
       ECC <- C0/lambda + C1*(-tau+n*T0+h*ARL2+d1*Tc+d2*Tr) + s*Cf/ARL1+Cr+(a+b*n)*(1/lambda-tau+n*T0+h*ARL2+d1*Tc+d2*Tr)/h
       cost <- ECC/ECT
     }else
     stop("You should at least give a pair of value to P0,P1 or C0,C1")
     return(cost)
   }
   
   mat <- array(NA,c(length(h),length(w),length(k)))
   for(i in 1:length(k))
     for(j in 1:length(w))
       mat[,j,i] = sapply(h,FUN=f,w=w[j],k=k[i],n=n)
   aa <- which(mat==min(mat),arr.ind=TRUE)
   optimum <- c(n,w[aa[1,2]],h[aa[1,1]],k[aa[1,3]],min(mat))
   names(optimum) <- c("n","Optimum w","Optimum h","Optimum k","Cost")
   
   par(mfrow=c(2,2),mar=c(5.1,4.1,1.1,1.1),oma=c(2,0,2,0))
   mat1 <- matrix(NA,length(h),length(k))
   for(j in 1:length(k))
     mat1[,j] <- sapply(h,FUN=f,w=optimum[2],k=k[j],n=optimum[1])
   contour(h,k,mat1,xlab="h",ylab="k",...)
   points(optimum[3],optimum[4],pch=3)
   
   mat2 <- matrix(NA,length(h),length(w))
   for(j in 1:length(w))
     mat2[,j] <- sapply(h,FUN=f,w=w[j],k=optimum[4],n=optimum[1])
   contour(h,w,mat2,xlab="h",ylab="w",...)
   points(optimum[3],optimum[2],pch=3)
   
   mat3 <- matrix(NA,length(k),length(w))
   for(j in 1:length(w))
     mat3[,j] <- sapply(k,FUN=f,w=w[j],h=optimum[3],n=optimum[1])
   contour(k,w,mat3,xlab="k",ylab="w",...)
   points(optimum[4],optimum[2],pch=3)
   mtext(sprintf('n=%s   Opt k=%s   Opt h=%s  Opt w=%s   Cost=%s',optimum[1], optimum[4], optimum[3],optimum[2],round(optimum[5],digits=4)),outer=T,line=.3,side=1)
   if(call.print){
     ca <- match.call()
     ca[[1]] <- NULL
     mtext(paste("ecoEwma(",paste(names(ca),"=",unlist(ca),sep="",collapse=", "),")"),side=3,line=0,outer=T)
   }
 }
