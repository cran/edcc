##' Calculate the optimum parameters, n(sample size), h(sampling
##' interval) and L(number of s.d. from control limits to center line)
##' for Econimic Design of X-bar control chart .
##'
##' For cost parameters P0, P1 and C0, C1, only one pair is needed.
##' If P0 and P1 are given, they will be used first, else C0 and C1
##' will be used.
##' For economic design of x-bar chart, only if the difference between
##' P0 and P1 keeps the same, the results are identical. If the
##' difference between C0 and C1 keeps the same, the optimum
##' parameters are identical but the cost values will change.
##' 
##' @title Economic design for X-bar chart
##' @param h sampling interval
##' @param L number of standard deviations from control limits to
##' center line.
##' @param n sample size.
##' @param lambda we assume the in-control time follows a exponential
##' distribution with mean 1/lambda. Default value is 0.05.
##' @param delta critical value: the extent of the shift when
##' assignable cause occurs. delta = |(mu1 - mu0)/sigma|. Default
##' value is 2.
##' @param P0 profit per hour earned by the process operating in
##' control. See 'Details'.
##' @param P1 profit per hour earned by the process operating out of
##' control
##' @param C0 cost per hour due to nonconformities produced while the
##' process is in control.
##' @param C1 cost per hour due to nonconformities produced while the
##' process is out of control.(C1 > C0)
##' @param Cr cost for searching and repairing the assignable cause,
##' including any downtime.
##' @param Cf cost per false alarm, including the cost of searching
##' for the cause and the cost of downtime if production ceases during
##' search.
##' @param T0 time to sample and chart one item.
##' @param Tc expected time to discover the assignable cause.
##' @param Tf expected search time when false alarm occures.
##' @param Tr expected time to repair the process.
##' @param a fixed cost per sample.
##' @param b cost per unit sampled.
##' @param d1 flag for whether production continues during searches
##' (1-yes, 0-no). Default value is 1.
##' @param d2 flag for whether production continues during repairs
##' (1-yes, 0-no). Default value is 1.
##' @param contour.plot a logical value indicating whether a contour
##' plot should be drawn. Default is FALSE.
##' @param call.print a logical value indicating whether the "call"
##' should be drawn on the contour plot. Default is TRUE
##' @param ... other arguments to be passed to contour function.
##' @return return the optimum parameters and the corresponding cost
##' value
##' @seealso \code{\link{ecoCusum}}, \code{\link{ecoEwma}}
##' @references Douglas (2009). Statistical quality control: a modern
##' introduction, sixth edition,  463-471.
##' 
##' Lorenzen and Vance (1986). The economic design of control charts:
##' a unified approach, Technometrics, 28. 3-10.
##' 
##' @examples 
##' # Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
##' # In control profit per hour is 110, out of control profit per hour is 10
##' x=ecoXbar(P0=110,P1=10,contour.plot=TRUE,nlevels=50)
##' summary(x)
##' contour(x)
##' # In control profit per hour is 150, out of control profit per hour
##' # is 50, the result is identical with the previous one, because the
##' #difference between P0 and P1 are the same
##' ecoXbar(P0=150,P1=50,contour.plot=TRUE,nlevels=50)
##' # In control cost per hour is 0, out of control cost per hour is 100.
##' #The result is the same with the previous one
##' ecoXbar(C0=0,C1=100,contour.plot=TRUE,nlevels=50)
##' # The optimum parameters are the same with the previous one,
##' # but Cost values are different. See 'details'
##' ecoXbar(C0=10,C1=110,contour.plot=TRUE,nlevels=50)
##' @export
ecoXbar <- function(h=seq(0.1,1,by=.01),L=seq(2,4.5,by=.01),n=1:15,lambda=.05,delta=2,P0=NULL,P1=NULL,C0=NULL,C1=NULL,Cr=25,Cf=50,T0=0.0167,Tc=1,Tf=0,Tr=0,a=1,b=.1,d1=1,d2=1,contour.plot=FALSE,call.print=TRUE,...){
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
  cost.frame <- NULL
  for(k in n){
    mat=outer(h,L,FUN=f,n=k)
    aa <- which(mat==min(mat),arr.ind=TRUE)
    cost.frame <- rbind(cost.frame,c(k,L[aa[1,2]],h[aa[1,1]],min(mat)))
  }
  colnames(cost.frame) <- c("n","Optimum L","Optimum h","Cost")
  rownames(cost.frame) <- rep("",length(n))
  optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
  if(contour.plot){
    par(mar=c(7.1,4.1,2.1,2.1))
    contour(h,L,outer(h,L,FUN=f,n=optimum[1]),xlab="h",ylab="L",...)
    points(optimum[3],optimum[2],pch=3)
    mtext(sprintf('n=%s   Opt L=%s   Opt h=%s   Cost=%s',optimum[1], optimum[2], optimum[3],round(optimum[4],digits=4)),side=1,line=4.5)
    if(call.print){
      ca <- match.call()
      ca[[1]] <- NULL
      mtext(paste("ecoXbar(",paste(names(ca),"=",unlist(ca),sep="",collapse=", "),")"),side=3,line=0)
  }
  }
  optXbar <- list(optimum=optimum,cost.frame=cost.frame)
  return(structure(optXbar,class="edcc",CALL=match.call()))
}


##' Calculate the optimum parameters of n(sample size), h(sampling
##' interval), k(reference value) and H(decision interval) for
##' Econimic Design of Cusum control chart. For more information about
##' the reference value see 'Details'.
##'
##' There is strong numerical and theoretical evidence that for given
##' L1, the value of L0 approaches its maximum when k(reference value)
##' is chosen mid-way the between AQL and the RQL: $k = mu0 +
##' 0.5*delta*sigma (Appl. Statist.(1974) 23, No. 3, p. 420). For this
##' reason we treat k as a constant value and optimize n, h and H.
##' For cost parameters P0, P1 and C0, C1, only one pair is needed.
##' If P0 and P1 are given, they will be used first, else C0 and C1
##' will be used.
##' 
##' @title Economic design for Cusum control chart
##' @param h sampling interval
##' @param H decision interval
##' @param n sample size
##' @param delta critical value: the extent of the shift when
##' assignable cause occurs. delta = |(mu1 - mu0)/sigma|. Default
##' value is 2.
##' @param lambda we assume the in-control time follows a exponential
##' distribution with mean 1/lambda. Default value is 0.05.
##' @param P0 profit per hour earned by the process operating in
##' control. See 'Details'.
##' @param P1 profit per hour earned by the process operating out of
##' control 
##' @param C0 cost per hour due to nonconformities produced while the
##' process is in control.
##' @param C1 cost per hour due to nonconformities produced while the
##' process is out of control.(C1 > C0)
##' @param Cr cost for searching and repairing the assignable cause,
##' including any downtime.
##' @param Cf cost per false alarm, including the cost of searching
##' for the cause and the cost of downtime if production ceases during
##' search.
##' @param T0 time to sample and chart one item.
##' @param Tc expected time to discover the assignable cause.
##' @param Tf expected search time when false alarm occures.
##' @param Tr expected time to repair the process.
##' @param a fixed cost per sample.
##' @param b cost per unit sampled.
##' @param d1 flag for whether production continues during searches
##' (1-yes, 0-no). Default value is 1.
##' @param d2 flag for whether production continues during repairs
##' (1-yes, 0-no). Default value is 1.
##' @param sided distinguish between one-, two-sided and Crosier's
##' modified two-sided CUSUM scheme by choosing "one", "two", and
##' "Crosier", respectively. See details in \code{\link[spc]{xcusum.arl}}
##' @param contour.plot a logical value indicating wether a contour
##' plot should be drawn. Default is FALSE.
##' @param call.print a logical value indicating whether the "call"
##' should be drawn on the contour plot. Default is TRUE
##' @param ... other arguments to be passed to contour function.
##' @return return the optimum parameters and the corresponding cost
##' value
##' @seealso \code{\link{ecoXbar}}, \code{\link{ecoEwma}},
##' \code{\link[spc]{xcusum.arl}}
##' @references Lorenzen and Vance (1986). The economic design of
##' control charts: a unified approach, Technometrics, 28. 3-10.
##' 
##' Taylor (1968). The economic design of cumulative sum control
##' charts, Technoinetrics, 10 479-488.
##' @examples
##' #Taylor (1968). Technoinetrics, 10, p427 Table3, row 1-4,14
##' y=ecoCusum(P0=150,P1=50,n=3:7,Cr=30,h=seq(1.3,1.5,by=0.01),d1=0,d2=0)
##' summary(y)
##' contour(y)
##' #ecoCusum(P0=150,P1=50,n=3:7,lambda=0.05,Cr=30,d1=0,d2=0)
##' #ecoCusum(P0=150,P1=140,n=3:7,h=seq(4,5,by=0.01),Cr=30,d1=0,d2=0)
##' #ecoCusum(P0=2000,P1=1000,n=4:8,Cr=30,h=seq(0.3,0.6,by=0.01),H=seq(0.8,1,by=0.01),d1=0,d2=0)
##' #ecoCusum(P0=150,P1=50,Cr=30,delta=0.5,h=seq(2.5,3,by=0.01),n=25:35,H=seq(0.2,0.4,by=0.01),d1=0,d2=0)
##' 
##' #Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
##' ecoCusum(h=seq(0.6,1,by=.01),H=seq(.1,1,by=.1),n=3:8,lambda=.05,
##' P0=110,P1=10,Cr=25,Cf=50,Tr=0,Tf=0,Tc=1,T0=.0167,a=1)
##' @export
ecoCusum <- function( h=seq(.1,2,by=.1), H=seq(.4,0.7,by=.01),n=3:7,delta = 2,lambda = .01, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 20, Cf = 10,T0 = 0, Tc = .1,Tf = .1,Tr = 0.2, a = .5, b = .1, d1 = 1, d2 = 1,sided = "one", contour.plot=FALSE,call.print=TRUE,...){
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
  cost.frame <- NULL
  for(nk in n){
    mat <- matrix(NA,length(h),length(H))
    for(i in 1:length(h))
      for(j in 1:length(H))
        mat[i,j] <- f(h[i],H[j],nk)
    aa <- which(mat==min(mat),arr.ind=TRUE)
    cost.frame <- rbind(cost.frame,c(nk,H[aa[1,][2]],h[aa[1,][1]],min(mat)))
  }
  colnames(cost.frame) <- c("n","Optimum H","Optimum h","Cost")
  rownames(cost.frame) <- rep("",length(n))
  optimum <- cost.frame[which(cost.frame[,4]==min(cost.frame[,4])),]
  if(contour.plot){
    par(mar=c(7.1,4.1,2.1,2.1))
    mat <- matrix(NA,length(h),length(H))
    for(i in 1:length(h))
      for(j in 1:length(H))
        mat[i,j] <- f(h[i],H[j],optimum[1])
    contour(h,H,mat,xlab="h",ylab="H",...)
    points(optimum[3],optimum[2],pch=3)
    mtext(sprintf('n=%s   Opt H=%s   Opt h=%s   Cost=%s',optimum[1], optimum[2], optimum[3],round(optimum[4],digits=4)),side=1,line=4.5)
    if(call.print){
      ca <- match.call()
      ca[[1]] <- NULL
      mtext(paste("ecoCusum(",paste(names(ca),"=",unlist(ca),sep="",collapse=", "),")"),side=3,line=0)
    }
  }
    optCusum <- list(optimum=optimum,cost.frame=cost.frame)
    return(structure(optCusum,class="edcc",CALL=match.call()))
}





##' Calculate the optimum parameters, n(sample size), h(sampling
##' interval), w(weight to the present sample) and k(number of
##' s.d. from control limits to center line) for econimic Design of
##' EWMA control chart .
##'
##' For cost parameters P0, P1 and C0, C1, only one pair is needed.
##' If P0 and P1 are given, they will be used first, else C0 and C1
##' will be used.
##' @title Economic design for EWMA control chart
##' @param h sampling interval
##' @param w the weight value given to the latest sample
##' @param k control limit coefficient
##' @param n sample size
##' @param delta critical value: the extent of the shift when
##' assignable cause occurs. delta = |(mu1 - mu0)/sigma|. Default
##' value is 2.
##' @param lambda we assume the in-control time follows a exponential
##' distribution with mean 1/lambda. Default value is 0.05.
##' @param P0 profit per hour earned by the process operating in
##' control. See 'Details'.
##' @param P1 profit per hour earned by the process operating out of
##' control 
##' @param C0 cost per hour due to nonconformities produced while the
##' process is in control.
##' @param C1 cost per hour due to nonconformities produced while the
##' process is out of control.(C1 > C0)
##' @param Cr cost for searching and repairing the assignable cause,
##' including any downtime.
##' @param Cf cost per false alarm, including the cost of searching
##' for the cause and the cost of downtime if production ceases during
##' search.
##' @param T0 time to sample and chart one item. 
##' @param Tc expected time to discover the assignable cause.
##' @param Tf expected search time when false alarm occures. 
##' @param Tr expected time to repair the process.
##' @param a fixed cost per sample.
##' @param b cost per unit sampled.
##' @param d1 flag for whether production continues during searches
##' (1-yes, 0-no). Default value is 1.
##' @param d2 flag for whether production continues during repairs
##' (1-yes, 0-no). Default value is 1.
##' @param sided distinguish between one- and two-sided EWMA control
##' chart by choosing "one" and "two", respectively. See details in
##' \code{\link[spc]{xewma.arl}}
##' @param contour.plot a logical value indicating whether a contour
##' plot should be drawn. Default is FALSE.
##' @param call.print a logical value indicating whether the "call"
##' should be drawn on the contour plot. Default is TRUE
##' @param ... other arguments to be passed to contour function.
##' @return return the optimum parameters and the corresponding cost
##' value
##' @seealso \code{\link{ecoXbar}}, \code{\link{ecoCusum}}
##' \code{\link[spc]{xewma.arl}}
##' @examples
##' #Douglas (2009). Statistical quality control: a modern introduction, sixth edition, p470.
##' x = ecoEwma(n=4:8,P0=110,P1=10,Cf=50,contour.plot=TRUE)
##' summary(x)
##' contour(x)
##' #ecoEwma(P0=150,P1=50,Cr=30,delta=0.5,h=seq(2.5,3,by=0.01),n=25:35)
##' @export
 ecoEwma <- function( h=seq(.7,1,by=.1), w=seq(0.7,1,by=.1),k=seq(2,4,by=0.1),n=4:8,delta = 2,lambda = .05, P0 = NULL, P1 = NULL,C0 = NULL,C1 = NULL, Cr = 25, Cf = 10,T0 = 0.0167,Tc = 1, Tf = 0, Tr = 0, a = 1, b = .1,d1=1,d2=1,sided="two",contour.plot=FALSE,call.print=TRUE,...){
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
   
   cost.frame <- NULL
   for(nk in n){
     mat <- array(NA,c(length(h),length(w),length(k)))
     for(i in 1:length(k))
       for(j in 1:length(w))
         mat[,j,i] = sapply(h,FUN=f,w=w[j],k=k[i],n=nk)
     aa <- which(mat==min(mat),arr.ind=TRUE)
     cost.frame <- rbind(cost.frame,c(nk,w[aa[1,2]],h[aa[1,1]],k[aa[1,3]],min(mat)))
   }
   colnames(cost.frame) <- c("n","Optimum w","Optimum h","Optimum k","Cost")
   rownames(cost.frame) <- rep("",length(n))
   optimum <- cost.frame[which(cost.frame[,5]==min(cost.frame[,5])),]
   if(contour.plot){
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
   optEwma <- list(optimum=optimum,cost.frame=cost.frame)
   return(structure(optEwma,class="edcc",CALL=match.call()))
 }


##' Print a "edcc" class object
##'
##' Print the optimum parameters of a "edcc" object
##' @title Print a "edcc" class object
##' @param x an object of class "edcc"
##' @param ... not used in this function
##' @return Invisible, the object itself.
##' @export
##' @S3method print edcc
##' @method print edcc
##' @examples x <- ecoXbar(P0=100,P1=0,nlevels=50)
##' print(x)
print.edcc <- function(x,...){
  print(x[[1]])
}

##' summary a "edcc" class object
##'
##' Both the optimum result and a data frame(cost.frame) are
##' returned. Each row of the cost.frame instead for an optimum result
##' corresponding to the specified n value.
##' @title summary a "edcc" class object
##' @param object an object of class "edcc"
##' @param ... not used in this function
##' @return A list containing the optimum reslut and a data frame is
##' returned.
##' @S3method summary edcc
##' @method summary edcc
##' @examples x <- ecoXbar(P0=100,P1=0,nlevels=50)
##' summary(x)
summary.edcc <- function(object,...){
  print(object[1]);print(object[2])
}


##' contour plot of "edcc" class
##'
##' contour plot for "edcc" class object
##' @title contour plot of "edcc" class
##' @param x an object of "edcc" class
##' @param ... arguments to be passed to contour plot, not supported
##' now
##' @return a contour plot
##' @S3method contour edcc
##' @method contour edcc
##' @examples z=ecoXbar(P0=100,P1=0,nlevels=50)
##' contour(z)
contour.edcc <- function(x,...){
  optimum <- x$optimum
  if(names(optimum)[2] == "Optimum L"){
    y <- attr(x,"CALL")
    y$n <- optimum[[1]]
    y$contour.plot <- NULL
    y[[1]] <- quote(.ecoxbar)
  }else
  if(names(optimum)[2] == "Optimum H"){
    y <- attr(x,"CALL")
    y$n <- optimum[[1]]
    y$contour.plot <- NULL
    y[[1]] <- quote(.ecocusum)
  }else
  if(names(optimum)[2] == "Optimum w"){
    y <- attr(x,"CALL")
    y$n <- optimum[[1]]
    y$contour.plot <- NULL
    y[[1]] <- quote(.ecoewma)
  }else
  cat("x should be an object returned by ecoXbar or ecoCusum or ecoEwma!\n")
  eval(y)
}

