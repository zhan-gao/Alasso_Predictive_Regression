# Prepare the plot
library(LasForecast)
library(lemon)
load("Result_obj_cv.RData")
H <- 3
m <- 4 # number of methods

df_trend <- NULL

change_bm <- function(coef_mat){
    
    method_names <- colnames(coef_mat)
    method_names[method_names == "b.m"] <- "bm"
    colnames(coef_mat) <- method_names
    
    return(coef_mat)
    
}

load("Monthly.koo.h0.0834.RData")
X <- change_bm(Monthly.Data[, -c(1,2,15)])
dates <- zoo::as.Date.yearmon(Monthly.Data$yyyymm)
X_plot <- reshape2::melt( cbind(dates, X) , id = "dates")

p_x_pred <- ggplot(data = X_plot) +
    geom_line(mapping = aes(x = dates,
                            y = value), color = "navyblue") +
    labs(x = NULL, y = NULL) +
    theme(
        panel.background =  element_blank(),
        panel.border = element_rect(
            linetype = 1,
            colour = "black",
            fill = NA
        ),
        panel.grid.major = element_line(linetype = 2, color = "grey90"),
        strip.background = element_blank(),
        strip.text =  element_text(face = "bold")
    ) +
    facet_wrap( ~ variable, ncol = 2, scales = "free")

p_x_pred_right <- p_x_pred + facet_wrap( ~ variable, ncol = 2, scales = "free", strip.position = "right")

ggsave(
    filename = "predictor_plot.png",
    plot = p_x_pred,
    device = "png",
    width = 170,
    height = 245,
    unit = "mm"
)

ggsave(
    filename = "predictor_plot_2.png",
    plot = p_x_pred_right,
    device = "png",
    width = 170,
    height = 245,
    unit = "mm"
)

for(h in c(0.0834, 0.25, 0.5, 1:H)){
    
    load(paste0("Monthly.koo.h", h, ".RData"))
    ifelse(h == 0.0834, h_print <- 0.083, h_print <- h)

    for (roll_window in c(120, 180)) {
        
        result_now <- eval(as.symbol(paste0("Result_", roll_window, "_", h, "_cv")))
        
        D <- Monthly.Data[!is.na(Monthly.Data$LongReturn), ]
        D <- change_bm(D)
        dates <- zoo::as.Date.yearmon(D$yyyymm[-(1:result_now$roll_window)])
        
        y_0 <- result_now$y[-(1:result_now$roll_window)]
        y_hat <- cbind(result_now$Lasso$y_hat,
                       result_now$Lasso_Std$y_hat,
                       result_now$ALasso$y_hat,
                       result_now$RepLasso$y_hat)
        colnames(y_hat) <- c("Plasso", "Slasso", "Alasso", "TAlasso")
        
        
        df_now <-
            cbind(data.frame(date = dates, y_0 = y_0),
                  as.data.frame(y_hat),
                  data.frame(h = paste0("h = ", h_print), rw = paste0(as.integer(roll_window/12), "-year Rolling Window")))
        df_now <- reshape2::melt(df_now, id = c("date", "h", "rw"))
        if(h %in% c(0.0834, 1, 3)){df_trend <- rbind(df_trend, df_now)}
        
        
        
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
            file = paste0("coef_plot_", roll_window, "_", h, ".png"),
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
                h_print
            )) + theme(plot.title = element_text(size = 10, hjust = 0.5))
        )
        
        ggsave(
            file = paste0("coef_plot_coint_dpdy", roll_window, "_", h, ".png"),
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
                h_print
            )) + theme(plot.title = element_text(size = 10, hjust = 0.5))
        )
        
        ggsave(
            file = paste0("coef_plot_coint_epdp", roll_window, "_", h, ".png"),
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
                h_print
            )) + theme(plot.title = element_text(size = 10, hjust = 0.5))
        )
        
        ggsave(
            file = paste0("coef_plot_coint_epdy", roll_window, "_", h, ".png"),
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

# -------------------------------------------------------------------------
# ----                          Trend plot                             ----
# -------------------------------------------------------------------------


methods_use <- c("Plasso", "Slasso", "Alasso", "TAlasso")
p_trend <- ggplot(data = df_trend) +
    geom_line(mapping = aes(
        x = date,
        y = value,
        color = variable,
        size = variable,
        alpha = variable
    )) +
    scale_colour_manual(breaks = c("y_0", methods_use),
                        values = c("gray", "red", "orange", "green", "blue"),
                        labels = c("True value", methods_use),
                        guide = guide_legend(override.aes = aes(alpha = NA))) +
    scale_size_manual(breaks = c("y_0", methods_use),
                      values =  rep(0.6, 5),
                      labels = c("True value", methods_use)) +
    scale_alpha_manual(breaks = c("y_0", methods_use),
                       values = c(0.5, 0.6, rep(0.5,3)),
                       labels = c("True value", methods_use)) +
    # scale_x_date(breaks = pt, labels = as.character(lubridate::year(pt))) +
    labs(x = NULL, y = NULL) +
    theme(
        panel.background =  element_blank(),
        panel.border = element_rect(
            linetype = 1,
            colour = "black",
            fill = NA
        ),
        panel.grid.major = element_line(linetype = 2, color = "grey90"),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text =  element_text(face = "bold")
    ) +
    facet_grid(vars(h), vars(rw), scales = "free")

ggsave(
    filename = "combined_trend.png",
    plot = p_trend,
    device = "png",
    width = 170,
    height = 245,
    unit = "mm"
)

ggsave(
    filename = "combined_trend_horizontal.png",
    plot = p_trend,
    device = "png",
    width = 245,
    height = 170,
    unit = "mm"
)

# -------------------------------------------------------------------------
# ----                     R2 bar plot - CV                            ----
# -------------------------------------------------------------------------

library(readxl)
library(ggplot2)

D <- read_excel("RFIT_read.xlsx")
D <- as.data.frame(t(D[, -(1:4)]))
method_names <- factor(
    c("Alasso", "TAlasso", "Plasso", "Slasso"),
    levels = c("Alasso", "TAlasso", "Plasso", "Slasso")
)

df_r2 <- NULL
case <- 0
for (roll_window in c(120, 180)){
    for(h in c(0.083, 0.25, 0.5, 1:3)){
        
        case <- case + 1
        df_now <- data.frame(x = method_names, y = D[, case])
        df_now <-
            cbind(df_now, data.frame(
                h = paste0("h = ", h),
                rw = paste0(as.integer(roll_window / 12), "-year Rolling Window")
            ))
        df_r2 <- rbind(df_r2, df_now)
        
    }
}

p_r2 <- ggplot(data = df_r2, aes(x = x, y = y, fill = x)) +
    geom_bar(stat="identity" ) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(ylim = c(0.6, 1)) +
    scale_fill_manual(
        values = c("Plasso" ="red", 
                   "Slasso" = "orange",
                   "Alasso" = "green", 
                   "TAlasso" = "blue")
    ) +
    theme(
        panel.background =  element_blank(),
        panel.border = element_rect(
            linetype = 1,
            colour = "black",
            fill = NA
        ),
        panel.grid.major = element_line(linetype = 2, color = "grey90"),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text =  element_text(face = "bold")
    ) +
    facet_grid(vars(h), vars(rw), scales = "free") + 
    guides(fill=FALSE)

ggsave(
    filename = "rfit_bar_60.png",
    plot = p_r2,
    device = "png",
    width = 170,
    height = 245,
    unit = "mm"
)

save.image(file = "plot_figures.RData")




