library(MASS)
library(robustbase)
library(matrixStats)


outlier_mahalanobis <- function(X, estimator = c("MCD", "MVE"), alpha = 0.975, tau = 1) {
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)
    estimator <- match.arg(estimator)

    covfit <- switch(estimator,
                     MCD = covmcd(X),
                     MVE = covmve(X)
    )
    center <- covfit$center
    covmat <- covfit$cov

    d <- sqrt(mahalanobis(X, center, covmat))
    u <- sqrt(qchisq(alpha, df = p))

    w1 <- ifelse(d <= u, 1, u / d)
    w2 <- w1^(2)/ tau

    mu_rew <- colSums(w1 * X) / sum(w1)

    Xc <- sweep(X, 2, mu_rew)
    cov_rew <- t(Xc) %*% (Xc * w2) / n

    outlier_flag <- d > u

    list(
        distance = d,
        # outlier = outlier_flag,
        out.id = which(outlier_flag),
        # keep = which(!outlier_flag),
        center = center,
        cov = covmat,
        # w1 = w1,
        # w2 = w2,
        mu_reweighted = mu_rew,
        cov_reweighted = cov_rew
    )
}

# random point
# 先跑这个更新outpro function它能return distance
outpro <- function(m,gval=NA,center=NA,plotit=FALSE,op=TRUE,MM=FALSE,cop=3,
                   xlab="VAR 1",ylab="VAR 2",STAND=TRUE,tr=FALSE,q=FALSE,pr=TRUE,...){
    m<-as.matrix(m)
    if(pr){
        if(!STAND){
            if(ncol(m)>1)print("STAND=FALSE. If measures are on different scales, might want to use STAND=TRUE")
        }}
    library(MASS)
    m=elimna(m)
    m<-as.matrix(m)
    nv=nrow(m)
    if(ncol(m)==1){
        dis<-(m-median(m,na.rm=TRUE))^2/mad(m,na.rm=TRUE)^2
        dis<-sqrt(dis)
        dis[is.na(dis)]=0
        crit<-sqrt(qchisq(.975,ncol(x)))
        chk<-ifelse(dis>crit,1,0)
        vec<-c(1:nrow(m))
        outid<-vec[chk==1]
        keep<-vec[chk==0]
    }
    if(ncol(m)>1){
        M=m
        if(STAND)m=standm(m,est=median,scat=mad)
        if(is.na(gval) && cop==1)gval<-sqrt(qchisq(.95,ncol(m)))
        if(is.na(gval) && cop!=1)gval<-sqrt(qchisq(.975,ncol(m)))
        if(cop==1 && is.na(center[1])){
            if(ncol(m)>2)center<-dmean(m,tr=.5,cop=1)
            if(ncol(m)==2){
                tempd<-NA
                for(i in 1:nrow(m))
                    tempd[i]<-depth(m[i,1],m[i,2],m)
                mdep<-median(tempd)
                flag<-(tempd==mdep)
                if(sum(flag)==1)center<-m[flag,]
                if(sum(flag)>1)center<-apply(m[flag,],2,mean)
            }}
        if(cop==2 && is.na(center[1])){
            center<-cov.mcd(m)$center
        }
        if(cop==4 && is.na(center[1])){
            center<-cov.mve(m)$center
        }
        if(cop==3 && is.na(center[1])){
            center<-apply(m,2,median)
        }
        if(cop==5 && is.na(center[1])){
            center<-tbs(m)$center
        }
        if(cop==6 && is.na(center[1])){
            center<-rmba(m)$center
        }
        if(cop==7 && is.na(center[1])){
            center<-spat(m)
        }
        flag<-rep(0, nrow(m))
        outid <- NA
        vec <- c(1:nrow(m))
        for (i in 1:nrow(m)){
            B<-m[i,]-center
            dis<-NA
            BB<-B^2
            bot<-sum(BB)
            if(bot!=0){
                for (j in 1:nrow(m)){
                    A<-m[j,]-center
                    temp<-sum(A*B)*B/bot
                    dis[j]<-sqrt(sum(temp^2))
                }
                temp<-idealf(dis)
                if(!MM)cu<-median(dis)+gval*(temp$qu-temp$ql)
                if(MM)cu<-median(dis)+gval*mad(dis)
                outid<-NA
                temp2<-(dis> cu)
                flag[temp2]<-1
            }}
        if(sum(flag) == 0) outid <- NA
        if(sum(flag) > 0)flag<-(flag==1)
        outid <- vec[flag]
        idv<-c(1:nrow(m))
        keep<-idv[!flag]
        if(ncol(m)==2){
            if(plotit){
                m=M # plot data using the original scale.
                plot(m[,1],m[,2],type="n",xlab=xlab,ylab=ylab)
                points(m[keep,1],m[keep,2],pch="*")
                if(length(outid)>0)points(m[outid,1],m[outid,2],pch="o")
                if(op){
                    tempd<-NA
                    keep<-keep[!is.na(keep)]
                    mm<-m[keep,]
                    for(i in 1:nrow(mm))tempd[i]<-depth(mm[i,1],mm[i,2],mm)
                    mdep<-median(tempd)
                    flag<-(tempd==mdep)
                    if(sum(flag)==1)center<-mm[flag,]
                    if(sum(flag)>1)center<-apply(mm[flag,],2,mean)
                    m<-mm
                }
                points(center[1],center[2],pch="+")
                x<-m
                temp<-fdepth(m,plotit=FALSE)
                flag<-(temp>=median(temp))
                xx<-x[flag,]
                xord<-order(xx[,1])
                xx<-xx[xord,]
                temp<-chull(xx)
                lines(xx[temp,])
                lines(xx[c(temp[1],temp[length(temp)]),])
            }}}
    list(n=nv,n.out=length(outid),out.id=outid,keep=keep,dis=dis, temp = temp, bot = bot)
}

## 然后用这个
# exhaustive
outlier_projection <- function(X, gval = NA, cop=2, MM = FALSE, alpha = 0.975, tau = 1, ...) {
    X <- as.matrix(X)
    n <- nrow(X)
    p <- ncol(X)

    outpro_result <- outpro(X, cop = cop, gval = gval, MM = MM, plotit = FALSE, ...)

    out.id <- outpro_result$out.id
    keep <- outpro_result$keep
    d <- outpro_result$dis  # vector of projection distances, length n

    u <- if (is.na(gval)) sqrt(qchisq(alpha, df = p)) else gval
    w1 <- ifelse(d <= u, 1, u / d)
    w2 <- w1^(2)/ tau

    # Reweighted mean
    mu_rew <- colSums(w1 * X) / sum(w1)

    # Reweighted covariance
    Xc <- sweep(X, 2, mu_rew)
    cov_rew <- t(Xc) %*% (Xc * w2) / n

    outlier_flag <- logical(n)
    if(length(out.id)) outlier_flag[out.id] <- TRUE

    list(
        outlier = outlier_flag,
        out.id = out.id,
        keep = keep,
        d = d,
        w1 = w1,
        w2 = w2,
        mu_reweighted = mu_rew,
        cov_reweighted = cov_rew
    )
}

