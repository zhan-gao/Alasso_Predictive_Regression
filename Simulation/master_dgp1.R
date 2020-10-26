# master file template

# -----------------------------------------------------
# Global environment setup for all DGP
# -----------------------------------------------------
rm(list = ls())

source_files <- c("dgp.R", "main_functions.R")
source_files <- lapply(source_files, source)
pacman::p_load(LasForecast, foreach, doParallel, doRNG)

t_range_temp <- c(80, 120, 200, 400, 800)
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
load("dgp_1_parameter.RData")
a0 <- dgp_1_parameter$a0
b0 <- dgp_1_parameter$b0
sigma_mat <- dgp_1_parameter$sigma_mat
phi <- dgp_1_parameter$phi
ur_index <- dgp_1_parameter$ur_index

# -----------------------------------------------------
# fix lambda
# -----------------------------------------------------
set.seed(100)
t0 <- Sys.time()
c.sel <- fix.lambda(
    100,
    100,
    a0,
    b0,
    generate_data_pure_unit_root,
    ur_index,
    c_seq = seq(0.1, 2, length.out = 1000),
    sigma_w_mat = sigma_mat,
    phi = phi,
    k = 5
)
t1 <- Sys.time() - t0

# -----------------------------------------------------
# run simulation
# -----------------------------------------------------
set.seed(99)
result <-
    run_simulation(
        generate_data_pure_unit_root,
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
cat("DGP 1: ", c.sel, "\n", file = "c_selection.txt", append = TRUE)
write.table(res_summary$fmse, file = "result_est_fmse.csv", append = TRUE, row.names = FALSE)
write.table(res_summary$coef, file = "result_est_coef.csv", append = TRUE, row.names = FALSE)
