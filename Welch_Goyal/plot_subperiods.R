# Plots by subperiods.

library(LasForecast)
library(lemon)
library(tidyverse)

load("Result_obj_cv.RData")
H <- 3
m <- 4 # number of methods

change_bm <- function(coef_mat){
    
    method_names <- colnames(coef_mat)
    method_names[method_names == "b.m"] <- "bm"
    colnames(coef_mat) <- method_names
    
    return(coef_mat)
    
}

# ---- Prepare dataframe -----
df_trend <- NULL

for(h in c(0.0834, 1, 3)){
    
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
        
        df_trend <- rbind(df_trend, df_now)
    }
}

df_trend_65 <- filter(df_trend, date > as.Date("1964-12-01") & date < as.Date("1975-01-01"))
df_trend_75 <- filter(df_trend, date > as.Date("1974-12-01") & date < as.Date("1985-01-01"))
df_trend_85 <- filter(df_trend, date > as.Date("1984-12-01") & date < as.Date("1995-01-01"))
df_trend_95 <- filter(df_trend, date > as.Date("1994-12-01") & date < as.Date("2005-01-01"))
df_trend_05 <- filter(df_trend, date > as.Date("2004-12-01"))

# ---- Plot ----

plot_trend_fn <- function(df, 
                          methods_use = c("Plasso", "Slasso", "Alasso", "TAlasso")){
    p_trend <- ggplot(data = df) +
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
    
    return(p_trend)
    
}

p_65 <- plot_trend_fn(df_trend_65)
p_75 <- plot_trend_fn(df_trend_75)
p_85 <- plot_trend_fn(df_trend_85)
p_95 <- plot_trend_fn(df_trend_95)
p_05 <- plot_trend_fn(df_trend_05)

save_plots <- function(f_name, f_name_h, p){
    
    ggsave(
        filename = f_name,
        plot = p,
        device = "png",
        width = 170,
        height = 245,
        unit = "mm"
    )
    
    ggsave(
        filename = f_name_h,
        plot = p,
        device = "png",
        width = 245,
        height = 170,
        unit = "mm"
    )
}

save_plots("sub_period_trend_65_75.png", 
           "sub_period_trend_65_75_horizontal.png", 
           p_65)
save_plots("sub_period_trend_75_85.png", 
           "sub_period_trend_75_85_horizontal.png", 
           p_75)
save_plots("sub_period_trend_85_95.png", 
           "sub_period_trend_685_95_horizontal.png", 
           p_85)
save_plots("sub_period_trend_95_05.png", 
           "sub_period_trend_95_05_horizontal.png", 
           p_95)
save_plots("sub_period_trend_05_12.png", 
           "sub_period_trend_05_12_horizontal.png", 
           p_05)
