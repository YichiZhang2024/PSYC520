data <- HolzingerSwineford1939[157:301, -c(1:6)]
rownames(data) <- 1:nrow(data)
# ## focus on one factor first
mod1 <- 'f1 =~ x1 + x2 + x3'
## add sex as a grouping variable
data_gp <- cbind(data, "sex" = HolzingerSwineford1939[157:301,]$sex)
## fit a one factor cfa model
fit_dat1 <- cfa(mod1, data = data_gp,
                group = "sex",
                estimator = "MLM",
                std.lv = TRUE)
## original data
res1 <- alignment(fit_dat1, group_name = c("female", "male"))
# CFA
path <- '
f1 =~ x1 + x2 + x3
f2 =~ x4 + x5 + x6
f3 =~ x7 + x8 + x9
'
model <- cfa(path, data = data, estimator = "MLM")
## builtin function from lavaan to extract lambda matrix
Lambda <- lavaan::inspect(model, what = "est")$lambda
## extract factor variance covariance matrix directly
Phi <- lavaan::inspect(model, what = "est")$psi

# Data set 2 -- Good leverage points
h = c(0.3881, 1.3762, 5.6153, 2.4312, 1.7442)
bartlett_predict <- lavPredict(model, method = "Bartlett")

data2_new_ob <- foreach(i = 141:145, j = 1:5) %do%{
    data[i, ] + h[[j]]*Lambda%*%bartlett_predict[i, ]
}

data2_new_ob <-data.frame(t(sapply(data2_new_ob,c)))

#replace the last 5 observations and the whole data set 2 is:
data2 <- rbind(data, data2_new_ob)
data2 <- data2[-(141:145), ]
rownames(data2) <- 1:nrow(data2)
data2 <- data.frame(matrix(unlist(data2), ncol=length(data2), byrow=FALSE))
data2 <- cbind(data2, "sex" = HolzingerSwineford1939[157:301,]$sex)
data2_fe <- data2 %>% filter(sex == 1)
data2_ma <- data2 %>% filter(sex == 2)

set.seed(3456)
res_mve <- robalign(method = "mve_mah", data_g1 = data2_fe[,1:3], data_g2 = data2_ma[,1:3], mod = mod2, group_name = c("female", "male"))
res_mcd <- robalign(method = "mcd_mah", data_g1 = data2_fe[,1:3], data_g2 = data2_ma[,1:3],
                    mod = mod2, group_name = c("Female", "Male"))
res_pro_mve <- robalign(method = "pro_mve", data_g1 = data2_fe[,1:3], data_g2 = data2_ma[,1:3], mod = mod2, group_name = c("female", "male"))
res_pro_mcd <- robalign(method = "pro_mcd", data_g1 = data2_fe[,1:3], data_g2 = data2_ma[,1:3], mod = mod2, group_name = c("female", "male"))
res_outmve <- robalign(method = "outpro_mve", data_g1 = data2_fe[,1:3], data_g2 = data2_ma[,1:3], mod = mod2, group_name = c("female", "male"))
res_outmcd <- robalign(method = "outpro_mcd", data_g1 = data2_fe[,1:3], data_g2 = data2_ma[,1:3], mod = mod2, group_name = c("female", "male"))

## covariance table
cov_mve <- as.data.frame(res_mve$rob_est_g1$weighted.covariance + res_mve$rob_est_g2$weighted.covariance)/2
cov_mcd <- as.data.frame(res_mcd$rob_est_g1$weighted.covariance + res_outmcd$rob_est_g2$weighted.covariance)/2
cov_outmve <- as.data.frame(res_outmve$rob_est_g1$weighted.covariance + res_outmve$rob_est_g2$weighted.covariance)/2
cov_outmcd <- as.data.frame(res_outmcd$rob_est_g1$weighted.covariance + res_outmcd$rob_est_g2$weighted.covariance)/2
cov_promve <- as.data.frame(res_pro_mve$rob_est_g1$weighted.covariance + res_pro_mve$rob_est_g2$weighted.covariance)/2
cov_promcd <- as.data.frame(res_pro_mcd$rob_est_g1$weighted.covariance + res_pro_mcd$rob_est_g2$weighted.covariance)/2
cov_tab <- cbind(cov(data_gp[1:3]), cov_mve, cov_mcd, cov_outmve, cov_outmcd, cov_promve, cov_promcd)
rownames(cov_tab) <- c("x1", "x2", "x3")
colnames(cov_tab) <- rep(paste0("x", 1:3), 7)

## factor loadings
ld_tab <- cbind(res1$lambda.aligned,
                res_mve$align_res$lambda.aligned, res_mcd$align_res$lambda.aligned,
                res_pro_mve$align_res$lambda.aligned, res_pro_mcd$align_res$lambda.aligned,
                res_outmve$align_res$lambda.aligned, res_outmcd$align_res$lambda.aligned)
rownames(ld_tab) <- c("Female", "Male")
colnames(ld_tab) <- rep(paste0("x", 1:3), 7)
## intercepts
int_tab <- cbind(res1$nu.aligned,
                 res_mve$align_res$nu.aligned, res_mcd$align_res$nu.aligned,
                 res_pro_mve$align_res$nu.aligned, res_pro_mcd$align_res$nu.aligned,
                 res_outmve$align_res$nu.aligned, res_outmcd$align_res$nu.aligned)
rownames(int_tab) <- c("Female", "Male")
colnames(int_tab) <- rep(paste0("x", 1:3), 7)
align_param <- list(`factor loadings` = ld_tab, `intercepts` = int_tab)

## factor means and factor variances
m_tab <- rbind(res1$pars[,1],
               res_mve$align_res$pars[,1],
               res_mcd$align_res$pars[,1],
               res_outmve$align_res$pars[,1],
               res_outmcd$align_res$pars[,1],
               res_pro_mve$align_res$pars[,1],
               res_pro_mcd$align_res$pars[,1])
var_tab <- rbind(res1$pars[,2],
                 res_mve$align_res$pars[,2],
                 res_mcd$align_res$pars[,2],
                 res_outmve$align_res$pars[,2],
                 res_outmcd$align_res$pars[,2],
                 res_pro_mve$align_res$pars[,2],
                 res_pro_mcd$align_res$pars[,2])
m_var_tab <- cbind(m_tab, var_tab)
colnames(m_var_tab) <- rep(c("Female", "Male"),2)
rownames(m_var_tab) <- c("Original Data",
                         "MD with MVE",
                         "MD with MCD",
                         "PMVE (random)",
                         "PMCD (random)",
                         "PMVE (exhaustive)",
                         "PMCD (exhaustive)")

