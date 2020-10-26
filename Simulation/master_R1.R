# master file template

# -----------------------------------------------------
# Global environment setup for all DGP
# -----------------------------------------------------
rm(list = ls())

source_files <- c("dgp.R", "main_functions.R")
source_files <- lapply(source_files, source)
pacman::p_load(LasForecast, foreach, doParallel)

t_range_temp <- c(40, 80, 120, 200, 400, 800)
t_range <- as.list(t_range_temp)
names(t_range) <- paste("t=", as.character(t_range_temp), sep = "")
R <- 10000
# m methods compared: ols, oracle ols, lasso, lasso_std, adalasso, replasso
m <- 6
# ml shrinkage methods compared: lasso, lasso_std, adalasso, replaso
ml <- 4
lasso_names <- c("adalasso", "replasso", "lasso", "lasso_std")
method_names <- c(
  "oracle",
  "ols",
  "adalasso",
  "replasso",
  "lasso",
  "lasso_std",
  paste("post_", lasso_names, sep = "")
)


# -----------------------------------------------------
# load parameters
# -----------------------------------------------------
load("dgp_2_parameter.RData")
a0 <- dgp_2_parameter$a0
b0 <- dgp_2_parameter$b0
b0[9] <- 0.1
sigma_mat <- dgp_2_parameter$sigma_mat
phi <- dgp_2_parameter$phi
ur_index <- dgp_2_parameter$ur_index
ur_index <- c(8)

# -----------------------------------------------------
# fix lambda
# -----------------------------------------------------
# set.seed(99)
# t0 <- Sys.time()
# c.sel <-
#     fix.lambda(
#         200,
#         100,
#         a0,
#         b0,
#         generate_data_mixed_roots,
#         ur_index,
#         c_seq = seq(1e-2, 1, length.out = 2000),
#         sigma_w_mat = sigma_mat,
#         phi = phi,
#         k = 10
#     )
# t1 <- Sys.time() - t0
# print(t1)
c.sel <- c(0.56492496, 0.09617309, 0.17863182, 0.30541521)


# -----------------------------------------------------
# run simulation
# -----------------------------------------------------
set.seed(99)
result <-
  run_simulation(
    generate_data_mixed_roots,
    sigma_mat,
    phi,
    a0,
    b0,
    t_range,
    R,
    m,
    ml,
    lasso_names,
    method_names,
    c.sel,
    ur_index
  )
res_summary <- sum_result(result)
coint_coef_summary <- sum_coint_coef(result, inactive_index = c(6, 7))
cat("DGP 2_R1: ", c.sel, "\n", file = "c_selection.txt", append = TRUE)
write.table(res_summary$fmse, file = "result_est_var_R1.csv", append = TRUE, row.names = FALSE)
write.table(res_summary$coef, file = "result_est_var_R1.csv", append = TRUE, row.names = FALSE)
write.table(coint_coef_summary, file = "result_est_var_R1.csv", append = TRUE, row.names = FALSE)

# Sel.est <- (Coef.est != 0)
# Sel.True <-
#   array(matrix(sel, R, p, byrow = TRUE), c(R, p, ml))
# Correct.ratio <- apply(Sel.True == Sel.est, 3, mean)
# Select.ratio <-
#   1 - colMeans(apply((Sel.True - Sel.est) == 1, c(1, 3), sum) / sum(sel))
# Screen.ratio <-
#   1 - colMeans(apply((Sel.True - Sel.est) == -1, c(1, 3), sum) / (p - sum(sel)))
# 
# coef.analysis <- list(correct = Correct.ratio,
#                       select = Select.ratio,
#                       screen = Screen.ratio)