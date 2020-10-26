library(LasForecast)
H = 3
MSE <- NULL

for(h in c(0.0834, 0.25, 0.5, 1:H)){

    load(paste0("Monthly.koo.h", h, ".RData"))
    
    for (Roll.Window in c(120, 180)) {
        
        D.temp = Monthly.Data[!is.na(Monthly.Data$LongReturn), ]
        
        XX = D.temp[, -c(1, 2, 15)]
        yy = D.temp$LongReturn
        
        result = LasForecast::roll_predict(XX,
                                           yy,
                                           Roll.Window,
                                           as.integer(h * 12),
                                           train_method_las = "cv")
        
        file_name <- paste0("Result_", as.character(Roll.Window), "_", as.character(h), "_bic")
        assign(file_name, result)
        
        MSE <- rbind(MSE, c(h, Roll.Window, result$mse))
    }
}

save.image(file = "Result_obj_cv.RData")
