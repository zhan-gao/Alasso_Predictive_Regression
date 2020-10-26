## Result summary
library(LasForecast)
library(tidyverse)
# -------------------------------------------------------------------
# ---------------------------- CV -----------------------------------
# -------------------------------------------------------------------

load("Result_obj_cv.RData")

MSE <- as.data.frame(MSE)
colnames(MSE)[1:2] <- c("h", "rw")
MSE <- MSE[order(MSE$rw, MSE$h), c(1, 2, 4:9)]

MAE <- NULL
freq_talasso_kill <- NULL
for(h in c(0.0834, 0.25, 0.5, 1:H)) {
    for (roll_window in c(120, 180)) {
        
        result_now <-
            eval(as.symbol(paste0(
                "Result_", roll_window, "_", h, "_cv"
            )))
        
        y_0 <- result_now$y[-(1:result_now$roll_window)]
        y_hat <- cbind(
            result_now$OLS$y_hat,
            result_now$RWwD$y_hat,
            result_now$ALasso$y_hat,
            result_now$RepLasso$y_hat,
            result_now$Lasso$y_hat,
            result_now$Lasso_Std$y_hat
        )
        colnames(y_hat) <- c("OLS", "RWwD", "Alasso", "TAlasso", "Plasso", "Slasso")
        
        MAE_now <- c(h = h, rw = roll_window, apply(y_hat, 2, function(x){mean(abs(x - y_0))}))
        MAE <- rbind(MAE, MAE_now)
        
        # Compare alasso and talasso in variable selection for cointegration pairs
        coef_alasso <- (result_now$ALasso$beta_hat[, c("dp", "dy", "ep", "dfr",  "tms")] != 0)
        coef_talasso <- (result_now$RepLasso$beta_hat[, c("dp", "dy", "ep", "dfr", "tms")] != 0)
        num_alasso_sel <- colSums(coef_alasso)
        
        freq_talasso_kill_now <- c(h = h, rw = roll_window, colSums((coef_alasso - coef_talasso) == 1) / num_alasso_sel)
        freq_talasso_kill <- rbind(freq_talasso_kill, freq_talasso_kill_now)
    }
}
MAE <- as.data.frame(MAE)
MAE <- MAE[order(MAE$rw, MAE$h),]

freq_talasso_kill <- as.data.frame(freq_talasso_kill)
freq_talasso_kill <- freq_talasso_kill[order(freq_talasso_kill$rw, freq_talasso_kill$h),]

write.csv(MSE, file = "MSE_cv.csv")
write.csv(MAE, file = "MAE_cv.csv")
write.csv(freq_talasso_kill, file = "freq_talasso_elim_cv.csv")

# --------------------------------------------------------------------
# ---------------------------- BIC -----------------------------------
# --------------------------------------------------------------------
load("Result_obj_bic.RData")

MSE_bic <- as.data.frame(MSE)
colnames(MSE_bic)[1:2] <- c("h", "rw")
MSE_bic <- filter(MSE_bic, h <= 3, rw <= 180) %>% arrange(rw)
MSE_bic <- MSE_bic[, c(1, 2, 4:9)]

MAE_bic <- NULL
freq_talasso_kill_bic <- NULL
H <- 3
for(h in c(0.0834, 0.25, 0.5, 1:H)) {
    for (roll_window in c(120, 180)) {
        
        result_now <-
            eval(as.symbol(paste0(
                "Result_", roll_window, "_", h, "_bic"
            )))
        
        y_0 <- result_now$y[-(1:result_now$roll_window)]
        y_hat <- cbind(
            result_now$OLS$y_hat,
            result_now$RWwD$y_hat,
            result_now$ALasso$y_hat,
            result_now$RepLasso$y_hat,
            result_now$Lasso$y_hat,
            result_now$Lasso_Std$y_hat
        )
        colnames(y_hat) <- c("OLS", "RWwD", "Alasso", "TAlasso", "Plasso", "Slasso")
        
        MAE_now <- c(h = h, rw = roll_window, apply(y_hat, 2, function(x){mean(abs(x - y_0))}))
        MAE_bic <- rbind(MAE_bic, MAE_now)
        
        # Compare alasso and talasso in variable selection for cointegration pairs
        coef_alasso <- (result_now$ALasso$beta_hat[, c("dp", "dy", "ep", "dfr",  "tms")] != 0)
        coef_talasso <- (result_now$RepLasso$beta_hat[, c("dp", "dy", "ep", "dfr", "tms")] != 0)
        num_alasso_sel <- colSums(coef_alasso)
        
        freq_talasso_kill_now <- c(h = h, rw = roll_window, colSums((coef_alasso - coef_talasso) == 1) / num_alasso_sel)
        freq_talasso_kill_bic <- rbind(freq_talasso_kill_bic, freq_talasso_kill_now)
        
    }
}

MAE_bic <- as.data.frame(MAE_bic)
MAE_bic <- MAE_bic[order(MAE_bic$rw, MAE_bic$h),]

freq_talasso_kill_bic <- as.data.frame(freq_talasso_kill_bic)
freq_talasso_kill_bic <- freq_talasso_kill_bic[order(freq_talasso_kill_bic$rw, freq_talasso_kill_bic$h),]

write.csv(MSE_bic, file = "MSE_bic.csv")
write.csv(MAE_bic, file = "MAE_bic.csv")
write.csv(freq_talasso_kill_bic, file = "freq_talasso_elim_bic.csv")