# source() #source Rand's functions
mve_mahalanobis <- function(x){
  library(MASS)
  x=as.matrix(x)
  p=ncol(x)
  n=nrow(x)
  if(ncol(x)>=2)
    mve<-cov.mve(x)
  if(ncol(x)==1){
    mve<-vector("list")
    mve$center<-median(x)
    mve$cov <- mad(x)^2
  }
  dep <- mahalanobis(x,mve$center,mve$cov)
  dist <- sqrt(dep)
  crit <- sqrt(qchisq(0.975,p))
  chklev <- ifelse(dist>crit,1,0)
  vec <- c(1:nrow(x))
  levid <- vec[chklev==1]

  varphi = length(levid)/nrow(x)
  prob = 1-varphi
  crit2 <- qchisq(0.975, p)
  tau = (p*pchisq(crit2,p+2)+ crit2*(1-prob))/p
  w1 <- c()
  w2 <- c()
  for (i in 1:nrow(x)){
    if (dist[i] <= crit){
      w1[i] <- 1.0
    }
    else
      w1[i] <- crit/dist[i]
  }
  for (i in 1:nrow(x)){
    if (dist[i] <= crit){
      w2[i] <- 1.0/tau
    }
    else
      w2[i] <- w1[i]^2/tau
  }

  mu_s = rep(0,p)
  reweighted_data <- as.data.frame(w1*as.matrix(x))
  # mu_s = colMeans(reweighted_data)
  mu_s = colSums(reweighted_data)/sum(w1)

  reweighted_data_sig <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
      diff_i <- matrix(x[i, ] - mu_s, ncol = 1)  # column vector
      reweighted_data_sig <- reweighted_data_sig + w2[i] * (diff_i %*% t(diff_i))
  }

  sigma_s <- reweighted_data_sig/n
  rownames(sigma_s) <- colnames(sigma_s) <- names(mu_s)

  list(levid = levid, levid.n = length(levid), weighted.mean = mu_s, weighted.covariance = sigma_s)
}

mcd_mahalanobis <- function(x){
  library(MASS)
  x=as.matrix(x)
  p=ncol(x)
  n=nrow(x)
  if(ncol(x)>=2)
    mcd<-cov.mcd(x)
  if(ncol(x)==1){
    mcd<-vector("list")
    mcd$center<-median(x)
    mcd$cov <- mad(x)^2
  }
  dep<-mahalanobis(x,mcd$center,mcd$cov)
  dist <- sqrt(dep)
  crit <- sqrt(qchisq(0.975,p))
  chklev <- ifelse(dist>crit,1,0)
  vec <- c(1:nrow(x))
  levid <- vec[chklev==1]

  varphi = length(levid)/nrow(x)
  prob = 1-varphi
  crit2 <- qchisq(0.975, p)
  tau = (p*pchisq(crit2,p+2)+ crit2*(1-prob))/p
  w1 <- c()
  w2 <- c()
  for (i in 1:nrow(x)){
    if (dist[i] <= crit){
      w1[i] <- 1.0
    }
    else
      w1[i] <- crit/dist[i]
  }
  for (i in 1:nrow(x)){
    if (dist[i] <= crit){
      w2[i] <- 1.0/tau
    }
    else
      w2[i] <- w1[i]^2/tau
  }

  mu_s = rep(0,p)
  reweighted_data <- as.data.frame(w1*as.matrix(x))
  mu_s = colSums(reweighted_data)/sum(w1)

  reweighted_data_sig <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
      diff_i <- matrix(x[i, ] - mu_s, ncol = 1)  # column vector
      reweighted_data_sig <- reweighted_data_sig + w2[i] * (diff_i %*% t(diff_i))
  }

  sigma_s <- reweighted_data_sig/n
  rownames(sigma_s) <- colnames(sigma_s) <- names(mu_s)

  list(levid = levid, levid.n = length(levid), weighted.mean = mu_s, weighted.covariance = sigma_s)
}



projection_mve <- function(x){
  x <- as.matrix(x)
  p = ncol(x)
  n = nrow(x)
  center <- cov.mve(x)$center
  A <-  as.data.frame(matrix(nrow = n, ncol = p))
  B <- as.data.frame(matrix(nrow = n, ncol = p))
  C <- as.data.frame(matrix(nrow = n^2, ncol = p))

  for (i in 1:nrow(x)){
    for (j in 1:nrow(x)){
    A[i, ] <- x[i, ] - center
    B[j, ] <- x[j, ] - center
    }
  }

  C = data.frame()
  for (i in 1:nrow(x)){
    for (j in 1:nrow(x)){
      results <- t((as.matrix(A[i, ])%*%t(as.matrix(B[j,]))/as.matrix(B[j, ])%*%t(as.matrix(B[j, ])))[1]*t(as.matrix(B[j, ])))
      df <- data.frame(results)
      C <- rbind(C, df)
    }
  }
  C <- as.matrix(C)
  D <- matrix(nrow = n^2, ncol = 1)
  for (i in 1:nrow(C)){
    D[i, ] <-sqrt(sum(C[i, ]^2))
  }
  D <- as.data.frame(D)

  D$fixed_i <- rep(c(1:n), each =n)
  D$pt <- rep(c(1:n), n)
  gval<-sqrt(qchisq(.975,ncol(x)))
  crit<-median(D$V1)+gval*(idealf(D$V1)$qu - idealf(D$V1)$ql)

  out <- D[which(D$V1 >= crit),]
  outlying_distances <- aggregate(out$V1, by = list(out$pt), median)
  levid <- unique(out$pt)

  varphi = length(levid)/nrow(x)
  prob = 1-varphi
  crit2 <- qchisq(0.975, p)
  tau = (p*pchisq(crit2,p+2)+ crit2*(1-prob))/p
  w1 <- as.data.frame(matrix(nrow=n, ncol=1))
  w2 <- as.data.frame(matrix(nrow=n, ncol=1))

  for (i in 1:nrow(w1)){
    w1[rownames(w1)[!(rownames(w1) %in% levid)],] <- 1
    for (j in 1:nrow(outlying_distances))
    if (rownames(w1)[i] == outlying_distances$Group.1[j]){
      w1[i,] <- crit/outlying_distances$x[j]
    }
  }

  for (i in 1:nrow(w2)){
    w2[rownames(w2)[!(rownames(w2) %in% levid)],] <- 1.0/tau
    for (j in 1:nrow(outlying_distances))
      if (rownames(w2)[i] == outlying_distances$Group.1[j]){
        w2[i,] <- as.matrix(w1)[i]^2/tau
      }
  }

  w1 <- as.numeric(as.matrix(w1))
  w2 <- as.numeric(as.matrix(w2))


  mu_s = rep(0,p)
  reweighted_data <- as.data.frame(w1*as.matrix(x))
  mu_s = colSums(reweighted_data)/sum(w1)

  reweighted_data_sig <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
      diff_i <- matrix(x[i, ] - mu_s, ncol = 1)  # column vector
      reweighted_data_sig <- reweighted_data_sig + w2[i] * (diff_i %*% t(diff_i))
  }

  sigma_s <- reweighted_data_sig/n
  rownames(sigma_s) <- colnames(sigma_s) <- names(mu_s)

  list(levid = levid, levid.n = length(levid), weighted.mean = mu_s, weighted.covariance = sigma_s)
}
#projection_mve(data2[,1:3])

projection_mcd <- function(x){
  x <- as.matrix(x)
  p = ncol(x)
  n = nrow(x)
  center <- cov.mcd(x)$center
  A <-  as.data.frame(matrix(nrow = n, ncol = p))
  B <- as.data.frame(matrix(nrow = n, ncol = p))
  C <- as.data.frame(matrix(nrow = n^2, ncol = p))

  for (i in 1:nrow(x)){
    for (j in 1:nrow(x)){
      A[i, ] <- x[i, ] - center
      B[j, ] <- x[j, ] - center
    }
  }

  C = data.frame()
  for (i in 1:nrow(x)){
    for (j in 1:nrow(x)){
      results <- t((as.matrix(A[i, ])%*%t(as.matrix(B[j,]))/as.matrix(B[j, ])%*%t(as.matrix(B[j, ])))[1]*t(as.matrix(B[j, ])))
      df <- data.frame(results)
      C <- rbind(C, df)
    }
  }
  C <- as.matrix(C)
  D <- matrix(nrow = n^2, ncol = 1)
  for (i in 1:nrow(C)){
    D[i, ] <-sqrt(sum(C[i, ]^2))
  }
  D <- as.data.frame(D)

  D$fixed_i <- rep(c(1:n), each =n)
  D$pt <- rep(c(1:n), n)
  gval<-sqrt(qchisq(.975,ncol(x)))
  crit<-median(D$V1)+gval*(idealf(D$V1)$qu - idealf(D$V1)$ql)

  out <- D[which(D$V1 >= crit),]
  outlying_distances <- aggregate(out$V1, by = list(out$pt), median)
  levid <- unique(out$pt)

  varphi = length(levid)/nrow(x)
  prob = 1-varphi
  crit2 <- qchisq(0.975, p)
  tau = (p*pchisq(crit2,p+2)+ crit2*(1-prob))/p
  w1 <- as.data.frame(matrix(nrow=n, ncol=1))
  w2 <- as.data.frame(matrix(nrow=n, ncol=1))

  for (i in 1:nrow(w1)){
    w1[rownames(w1)[!(rownames(w1) %in% levid)],] <- 1
    for (j in 1:nrow(outlying_distances))
      if (rownames(w1)[i] == outlying_distances$Group.1[j]){
        w1[i,] <- crit/outlying_distances$x[j]
      }
  }

  for (i in 1:nrow(w2)){
    w2[rownames(w2)[!(rownames(w2) %in% levid)],] <- 1.0/tau
    for (j in 1:nrow(outlying_distances))
      if (rownames(w2)[i] == outlying_distances$Group.1[j]){
        w2[i,] <- as.matrix(w1)[i]^2/tau
      }
  }

  w1 <- as.numeric(as.matrix(w1))
  w2 <- as.numeric(as.matrix(w2))


  mu_s = rep(0,p)
  reweighted_data <- as.data.frame(w1*as.matrix(x))
  mu_s = colSums(reweighted_data)/sum(w1)

  reweighted_data_sig <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
      diff_i <- matrix(x[i, ] - mu_s, ncol = 1)  # column vector
      reweighted_data_sig <- reweighted_data_sig + w2[i] * (diff_i %*% t(diff_i))
  }

  sigma_s <- reweighted_data_sig/n
  rownames(sigma_s) <- colnames(sigma_s) <- names(mu_s)

  list(levid = levid, levid.n = length(levid), weighted.mean = mu_s, weighted.covariance = sigma_s)
}

outpro_mve <- function(x){
  x=as.matrix(x)
  p=ncol(x)
  n=nrow(x)
  levid <- outpro_dis(x, cop=4)$out.id

  varphi = length(levid)/nrow(x)
  prob = 1-varphi
  crit <- qchisq(0.975, p)
  tau = (p*pchisq(crit,p+2)+ crit*(1-prob))/p
  # crit2 <- qchisq(0.975, p)
  # tau = (p*pchisq(crit2,p+2)+ crit2*(1-prob))/p
  w1 <- c()
  w2 <- c()
  dist <- outpro_dis(x, cop=4)$dis
  for (i in 1:nrow(x)){
    if (dist[i] <= crit){
      w1[i] <- 1.0
    }
    else
      w1[i] <- crit/dist[i]
  }
  for (i in 1:nrow(x)){
    if (dist[i] <= crit){
      w2[i] <- 1.0/tau
    }
    else
      w2[i] <- w1[i]^2/tau
  }

  mu_s = rep(0,p)
  reweighted_data <- as.data.frame(w1*as.matrix(x))
  mu_s = colSums(reweighted_data)/sum(w1)

  reweighted_data_sig <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
      diff_i <- matrix(x[i, ] - mu_s, ncol = 1)  # column vector
      reweighted_data_sig <- reweighted_data_sig + w2[i] * (diff_i %*% t(diff_i))
  }

  sigma_s <- reweighted_data_sig/n
  rownames(sigma_s) <- colnames(sigma_s) <- names(mu_s)

  list(levid = levid, levid.n = length(levid), weighted.mean = mu_s, weighted.covariance = sigma_s)
}

outpro_mcd <- function(x){
  x=as.matrix(x)
  p=ncol(x)
  n=nrow(x)
  levid <- outpro_dis(x, cop=2)$out.id

  varphi = length(levid)/nrow(x)
  prob = 1-varphi
  # crit2 <- qchisq(0.975, p)
  # tau = (p*pchisq(crit2,p+2)+ crit2*(1-prob))/p
  crit <- qchisq(0.975, p)
  tau = (p*pchisq(crit,p+2)+ crit*(1-prob))/p
  w1 <- c()
  w2 <- c()
  dist <- outpro_dis(x, cop=2)$dis
  for (i in 1:nrow(x)){
    if (dist[i] <= crit){
      w1[i] <- 1.0
    }
    else
      w1[i] <- crit/dist[i]
  }
  for (i in 1:nrow(x)){
    if (dist[i] <= crit){
      w2[i] <- 1.0/tau
    }
    else
      w2[i] <- w1[i]^2/tau
  }

  mu_s = rep(0,p)
  reweighted_data <- as.data.frame(w1*as.matrix(x))
  mu_s = colSums(reweighted_data)/sum(w1)

  reweighted_data_sig <- matrix(0, nrow = p, ncol = p)
  for (i in 1:n) {
      diff_i <- matrix(x[i, ] - mu_s, ncol = 1)  # column vector
      reweighted_data_sig <- reweighted_data_sig + w2[i] * (diff_i %*% t(diff_i))
  }

  sigma_s <- reweighted_data_sig/n
  rownames(sigma_s) <- colnames(sigma_s) <- names(mu_s)

  list(levid = levid, levid.n = length(levid), weighted.mean = mu_s, weighted.covariance = sigma_s)
}

## helper function extract parameter estimates
ext_est <- function(mod, rob_est, data) {
    fit <- cfa(
        mod,
        sample.cov  = rob_est$cov_reweighted,
        sample.mean = rob_est$mu_reweighted,
        sample.nobs = nrow(data),
        std.lv      = TRUE
    )
    list(
        lambda = t(lavInspect(fit, "est")$lambda),
        nu     = t(lavInspect(fit, "est")$nu)
    )
}
## robust alignment function
robalign <- function(method, data, mod, group_names){
    grouped_data <- split(data[, !names(data) %in% "group"], data$group)
    grouped_data <- grouped_data[group_names]
    if(method == "mve_mah"){
        rob_estimates <- lapply(grouped_data, function(grouped_data){
            outlier_mahalanobis(grouped_data[,3:5], estimator = "MVE")
        })
    }else if(method == "mcd_mah"){
        rob_estimates <- lapply(grouped_data, function(grouped_data){
            outlier_mahalanobis(grouped_data[,3:5], estimator = "MCD")
        })
    }else if(method == "pro_mve"){
        rob_estimates <- lapply(grouped_data, function(grouped_data){
            outlier_projection(X = grouped_data[,3:5], cop = 4)
        })
    }else if(method == "pro_mcd"){
        rob_estimates <- lapply(grouped_data, function(grouped_data){
            outlier_projection(X = grouped_data[,3:5], cop = 2)
        })
    # }else if(method == "outpro_mve"){
    #     rob_estimates <- lapply(grouped_data, function(grouped_data){
    #         outpro_mcd(x = grouped_data[,3:5], cov = 4)
    #     })
    # }else if(method == "outpro_mcd"){
    #     rob_estimates <- lapply(grouped_data, function(grouped_data){
    #         outpro_mcd(x = grouped_data[,3:5], cov = 2)
    #     })
    }

    fits <- list(
        ext_est(mod, rob_estimates[[1]], grouped_data[[1]]),
        ext_est(mod, rob_estimates[[2]], grouped_data[[2]]),
        ext_est(mod, rob_estimates[[3]], grouped_data[[3]]),
        ext_est(mod, rob_estimates[[4]], grouped_data[[4]])

    )
    # Combine results
    ld <- do.call(rbind, lapply(fits, `[[`, "lambda"))
    nu <- do.call(rbind, lapply(fits, `[[`, "nu"))

    # Assign row names
    rownames(ld) <- rownames(nu) <- names(grouped_data)

    res <- invariance.alignment(
        lambda = ld,
        nu = nu,
        wgt = matrix(sqrt(sapply(grouped_data, nrow)), nrow = 4, ncol = 3)
    )
    list(align_res = res, rob_est= rob_estimates)
}
## test
# group_names <- c("female_Pasteur", "Male_Pasteur", "female_Grant-White", "Male_Grant-White")
# res <- robalign(method = "mve_mah", data = data2, mod, group_names)
# res2 <- robalign(method = "pro_mve", data = data2, mod, group_names)

alignment <- function(model_fit, group_name){
    ld_dat <- rbind(t(lavInspect(model_fit,  what = "est")$`1`$lambda),
                    t(lavInspect(model_fit,  what = "est")$`2`$lambda))
    rownames(ld_dat) <- group_name
    ## reformat intercepts matrix
    int_dat <- rbind(t(lavInspect(model_fit,  what = "est")$`1`$nu),
                     t(lavInspect(model_fit,  what = "est")$`2`$nu))
    rownames(int_dat) <- group_name
    ## run the alignment function from sirt
    res <- invariance.alignment(
        lambda = ld_dat,
        nu = int_dat,
        wgt = matrix(sqrt(summary(model_fit)[[3]]$nobs), nrow = 2, ncol = 3)
    )
    res
}
# # test

#
outpro_dis <- function(m,gval=NA,center=NA,plotit=FALSE,op=TRUE,MM=FALSE,cop=3,
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

# halfdis <- function(x){
#   x=as.matrix(x)
#   #x <- standm(x)
#   p=ncol(x)
#   n=nrow(x)
#
#   dist<-prodepth(x)
#   crit <- sqrt(qchisq(0.975,p)) # FIX THIS
#   chklev <- ifelse(dist>crit,1,0)
#   vec <- c(1:nrow(x))
#   levid <- vec[chklev==1]
#
#   varphi = length(levid)/nrow(x)
#   prob = 1-varphi
#   crit2 <- qchisq(0.975, p)
#   tau = (p*pchisq(crit2,p+2)+ crit2*(1-prob))/p
#   w1 <- c()
#   w2 <- c()
#   for (i in 1:nrow(x)){
#     if (dist[i] <= crit){
#       w1[i] <- 1.0
#     }
#     else
#       w1[i] <- crit/dist[i]
#   }
#   for (i in 1:nrow(x)){
#     if (dist[i] <= crit){
#       w2[i] <- 1.0/tau
#     }
#     else
#       w2[i] <- w1[i]^2/tau
#   }
#
#   mu_s = rep(0,p)
#   reweighted_data <- as.data.frame(w1*as.matrix(data))
#   mu_s = colMeans(reweighted_data)
#
#   diff <- as.data.frame(matrix(nrow = n, ncol = 9, byrow=TRUE))
#   for (i in 1:n) {
#     diff[i, ] <- data[i, ] - mu_s
#   }
#   reweighted_data_sig <- as.data.frame(w2*as.matrix(diff))
#   sigma_s <- cov(reweighted_data_sig)
#
#   list(levid = levid, levid.n = length(levid), weighted.mean = mu_s, weighted.covariance = sigma_s)
# }



# projection_method <- function(x){
#   library(MASS)
#   x=as.matrix(x)
#   n<-nrow(x)
#   p=ncol(x)
#   dist<-1/prodepth(x)
#   crit <- median(dist) + sqrt(qchisq(0.975,p))*(idealf(dist)$qu - idealf(dist)$ql)
#   chklev <- ifelse(dist>crit,1,0)
#   vec <- c(1:nrow(x))
#   levid <- vec[chklev==1]
#
#   varphi = length(levid)/nrow(x)
#   prob = 1-varphi
#   crit2 <- qchisq(0.975, p)
#   tau = (p*pchisq(crit2,p+2)+ crit2*(1-prob))/p
#   w1 <- c()
#   w2 <- c()
#   for (i in 1:nrow(x)){
#     if (dist[i] <= crit){
#       w1[i] <- 1.0
#     }
#     else
#       w1[i] <- crit/dist[i]
#   }
#   for (i in 1:nrow(x)){
#     if (dist[i] <= crit){
#       w2[i] <- 1.0/tau
#     }
#     else
#       w2[i] <- w1[i]^2/tau
#   } ##到这里w1和w2都是对的 ##
#
#   mu_s = rep(0,p)
#   reweighted_data <- as.data.frame(w1*as.matrix(data))
#   mu_s = colMeans(reweighted_data)
#
#   diff <- as.data.frame(matrix(nrow = n, ncol = 9, byrow=TRUE))
#   for (i in 1:n) {
#     diff[i, ] <- data[i, ] - mu_s
#   }
#   reweighted_data_sig <- as.data.frame(w2*as.matrix(diff))
#   sigma_s <- cov(reweighted_data_sig)
#
#   list(levid = levid, levid.n = length(levid), weighted.mean = mu_s, weighted.covariance = sigma_s)
# }
