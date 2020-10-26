rm(list = ls())
load("GW_data.RData")
x <- dplyr::select(Monthly.Data, -c("yyyymm", "equity_premium", "LongReturn"))
y <- dplyr::select(Monthly.Data, "equity_premium")
names(x)[1] <- "bm"

x_mat <- as.matrix(x)
y_mat <- as.matrix(y)
u <- y_mat - (x_mat %*% solve(t(x_mat) %*% x_mat) %*% t(x_mat)) %*% y_mat
u <- u[-c(1,2)]

n <- length(u)

ar_resid <- sapply(x, function(x){
    ar_obj <- ar(x, aic = FALSE, order.max = 1, method = "yule-walker")
    ar_obj$resid[-1]
})

ar_coef <- apply(x, 2, function(x){
    ar_obj <- ar(x, aic = FALSE, order.max = 1, method = "yule-walker")$ar
})
c_lur <- (1 - ar_coef[!(names(ar_coef) %in% c("ltr", "infl", "svar"))]) * n

z_1 <- ar_resid[-1, "ltr"]
z_2 <- ar_resid[-1, "infl"]
z_3 <- ar_resid[-1, "svar"]

v_1 <- ar_resid[, 'dy']
v_2 <- lm(0 + x$dy ~ x$dp)$residuals[-1]
v_3 <- ar_resid[, 'tms']
v_4 <- lm(0 + x$tms ~ x$dfr)$residuals[-1]
ar_v <- apply(cbind(v_1, v_2, v_3, v_4), 2, function(x){
    ar(x, order.max = 1, aic = FALSE, method = "yule-walker")$ar
})

v_1_resid <- ar(v_1, aic = FALSE, order.max = 1, method = "yule-walker")$resid[-1]
v_2_resid <- ar(v_2, aic = FALSE, order.max = 1, method = "yule-walker")$resid[-1]
v_3_resid <- ar(v_3, aic = FALSE, order.max = 1, method = "yule-walker")$resid[-1]
v_4_resid <- ar(v_4, aic = FALSE, order.max = 1, method = "yule-walker")$resid[-1]

e_1 <- ar_resid[-1, "ep"]
e_2 <- ar_resid[-1, "bm"]
e_3 <- ar_resid[-1, "tbl"]
e_4 <- ar_resid[-1, "dfy"]
e_5 <- ar_resid[-1, "ntis"]
e_6 <- ar_resid[-1, "dy"]
e_7 <- ar_resid[-1, 'dp']
e_8 <- ar_resid[-1, "dfr"]
e_9 <- ar_resid[-1, "tms"]

# -----------------------------------------------------
# DGP 1 pure unit root
# -----------------------------------------------------
innovation_mat_1 <- cbind(e_1, e_2, e_3, e_4, e_5, e_6, e_7, e_8, e_9, u)
cov_temp <- cov(innovation_mat_1)

cov_mat_1 <- diag(diag(cov_temp)^{-1/2}) %*% cov_temp %*% diag(diag(cov_temp)^{-1/2})

dgp_1_parameter <- list(a0 = 0.25,
                        b0 = c(rep(1, 3), rep(0, 6)),
                        sigma_mat = cov_mat_1,
                        ur_index = 1:3)
save(dgp_1_parameter, file = "dgp_1_parameter.RData")

# -----------------------------------------------------
# DGP 2 mixed roots
# -----------------------------------------------------
innovation_mat_2 <- cbind(z_1, z_2, z_3, v_1_resid, v_2_resid, v_3_resid, v_4_resid, e_1, e_2, e_3, e_4, e_5, u)

cov_temp <- cov(innovation_mat_2)
cov_temp[nrow(cov_temp), 1:7] <- 0
cov_temp[1:7, ncol(cov_temp)] <- 0
cov_mat_2 <- diag(diag(cov_temp)^{-1/2}) %*% cov_temp %*% diag(diag(cov_temp)^{-1/2})

phi_z_2 <- ar_coef[c("ltr", "infl", "svar")]
names(phi_z_2) <- paste0("z_", 1:3)
phi_v_2 <- ar_v
phi_x_2 <- rep(0, 5)
names(phi_x_2) <- paste0("e_", 1:5)
phi_u_2 <- 0
names(phi_u_2) <- "u"
phi_2 <- c(phi_z_2, phi_v_2, phi_x_2, phi_u_2)

dgp_2_parameter <- list(a0 = 0.3,
                        b0 = c(0.4, 0, 0, 0.5, -0.5, 0, 0, 1, 1, 0, 0, 0),
                        sigma_mat = cov_mat_2,
                        phi = phi_2,
                        ur_index = c(8, 9))
save(dgp_2_parameter, file = "dgp_2_parameter.RData")

# -----------------------------------------------------
# DGP 3 mixed roots with local unit root (lur)
# -----------------------------------------------------
# ep(1), bm(2), tbl(3): local unit root
# ntis(4) and dfy(5): pure unit-root
dgp_3_parameter <- dgp_2_parameter

dgp_3_parameter$b0 <- c(0.4, 0, 0, 0.5, -0.5, 0, 0, 1, 0, 0, 1, 0)
dgp_3_parameter$ur_index <- c(8, 11)
save(dgp_3_parameter, file = "dgp_3_parameter.RData")

# -----------------------------------------------------
# DGP 4 ARDL
# -----------------------------------------------------
dgp_4_parameter <- dgp_2_parameter
dgp_4_b0 <- c(0.4, 0.75, -0.75, rep(0, 6), 1.5, 1, 0, 1, 0, 0, -1, 0, 0, 0, 0.6, 0.8, 0, 0.4, 0, 0)

dgp_4_parameter$b0 <- dgp_4_b0
dgp_4_parameter$ur_index <- c(10, 13)
save(dgp_4_parameter, file = "dgp_4_parameter.RData")