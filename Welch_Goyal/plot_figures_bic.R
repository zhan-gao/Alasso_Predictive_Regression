# Prepare the plot
library(LasForecast)
library(lemon)
load("Result_obj_bic.RData")
H <- 3
m <- 4 # number of methods

change_bm <- function(coef_mat){
    
    method_names <- colnames(coef_mat)
    method_names[method_names == "b.m"] <- "bm"
    colnames(coef_mat) <- method_names
    
    return(coef_mat)
    
}

for(h in c(0.0834, 0.25, 0.5, 1:H)){
    
    load(paste0("Monthly.koo.h", h, ".RData"))
    ifelse(h == 0.0834, h_print <- 0.083, h_print <- h)
    
    for (roll_window in c(120, 180)) {
        
        result_now <- eval(as.symbol(paste0("Result_", roll_window, "_", h, "_bic")))
        
        D <- Monthly.Data[!is.na(Monthly.Data$LongReturn), ]
        D <- change_bm(D)
        dates <- zoo::as.Date.yearmon(D$yyyymm[-(1:result_now$roll_window)])
        
        # -------- Coefficient plot ----------
        
        coef_plasso <- change_bm(result_now$Lasso$beta_hat)
        coef_slasso <- change_bm(result_now$Lasso_Std$beta_hat)
        coef_alasso <- change_bm(result_now$ALasso$beta_hat)
        coef_talasso <- change_bm(result_now$RepLasso$beta_hat)
        
        coef_est <- list(Plasso =  coef_plasso,
                         Slasso = coef_slasso,
                         Alasso = coef_alasso,
                         TAlasso = coef_talasso)
        
        assign(
            paste0("coef_plot_", roll_window, "_h_", h),
            plot_coef(coef_est, dates,
                      col_vec = c("green","red", "orange", "blue"),
                      line_size = rep(0.7, 4),
                      alpha_size = rep(0.55, 4),
                      num_col = 2)
        )
        
        ggsave(
            file = paste0("coef_plot_", roll_window, "_", h, "_bic.png"),
            plot = eval(as.symbol(
                paste0("coef_plot_", roll_window, "_h_", h)
            )),
            device = "png",
            width = 170,
            height = 245,
            unit = "mm"
        )
        
        # -------- Cointegrate pair plot ----------
        
        # dy - dp
        coef_est_coint_dpdy <-
            list(Alasso = coef_alasso[, c("dp", "dfr", "dy",  "tms")],
                 TAlasso = coef_talasso[, c("dp", "dfr", "dy",  "tms")])
        
        assign(
            paste0("coef_plot_coint_dpdy", roll_window, "_h_", h),
            plot_coef(
                coef_est_coint_dpdy,
                dates,
                col_vec = c("green", "blue"),
                line_size = rep(1, 2),
                alpha_size = rep(0.55, 2),
                num_col = 2
            ) + ggtitle(paste0(
                as.integer(roll_window / 12),
                "-year Rolling Window; ",
                "Horizon = ",
                h
            )) + theme(plot.title = element_text(size = 10, hjust = 0.5))
        )
        
        ggsave(
            file = paste0("coef_plot_coint_dpdy", roll_window, "_", h, "_bic.png"),
            plot = eval(as.symbol(
                paste0("coef_plot_coint_dpdy", roll_window, "_h_", h)
            )),
            device = "png",
            width = 245,
            height = 170,
            unit = "mm"
        )
        
        # ep - dp
        coef_est_coint_epdp <-
            list(Alasso = coef_alasso[, c("ep", "dfr", "dp",  "tms")],
                 TAlasso = coef_talasso[, c("ep", "dfr", "dp",  "tms")])
        
        assign(
            paste0("coef_plot_coint_epdp", roll_window, "_h_", h),
            plot_coef(
                coef_est_coint_epdp,
                dates,
                col_vec = c("green", "blue"),
                line_size = rep(1, 2),
                alpha_size = rep(0.55, 2),
                num_col = 2
            ) + ggtitle(paste0(
                as.integer(roll_window / 12),
                "-year Rolling Window; ",
                "Horizon = ",
                h
            )) + theme(plot.title = element_text(size = 10, hjust = 0.5))
        )
        
        ggsave(
            file = paste0("coef_plot_coint_epdp", roll_window, "_", h, "_bic.png"),
            plot = eval(as.symbol(
                paste0("coef_plot_coint_epdp", roll_window, "_h_", h)
            )),
            device = "png",
            width = 245,
            height = 170,
            unit = "mm"
        )
        
        # ep - dy
        coef_est_coint_epdy <-
            list(Alasso = coef_alasso[, c("ep", "dfr", "dy",  "tms")],
                 TAlasso = coef_talasso[, c("ep", "dfr", "dy",  "tms")])
        
        assign(
            paste0("coef_plot_coint_epdy", roll_window, "_h_", h),
            plot_coef(
                coef_est_coint_epdy,
                dates,
                col_vec = c("green", "blue"),
                line_size = rep(1, 2),
                alpha_size = rep(0.55, 2),
                num_col = 2
            ) + ggtitle(paste0(
                as.integer(roll_window / 12),
                "-year Rolling Window; ",
                "Horizon = ",
                h
            )) + theme(plot.title = element_text(size = 10, hjust = 0.5))
        )
        
        ggsave(
            file = paste0("coef_plot_coint_epdy", roll_window, "_", h, "_bic.png"),
            plot = eval(as.symbol(
                paste0("coef_plot_coint_epdy", roll_window, "_h_", h)
            )),
            device = "png",
            width = 245,
            height = 170,
            unit = "mm"
        )
        
    }
}