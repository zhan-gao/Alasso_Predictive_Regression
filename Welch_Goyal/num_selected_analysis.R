rm(list = ls())
load("Result_obj_cv.RData")

result_ratio <- NULL
for (rw in c(120, 180)) {
    for (h in c(0.0834, 0.25, 0.5, 1, 2, 3)) {
        variable_name <- paste0("Result_", rw, "_", h, "_cv")
        assign("result_temp", eval(as.symbol(variable_name)))
        result_temp$Lasso$beta_hat
        result_ratio <- rbind(result_ratio, 
                              c(mean(rowSums(result_temp$ALasso$beta_hat != 0) > 1),
                                mean(rowSums(result_temp$RepLasso$beta_hat != 0) > 1),
                                mean(rowSums(result_temp$Lasso$beta_hat != 0) > 1),
                                mean(rowSums(result_temp$Lasso_Std$beta_hat != 0) > 1)
                              ))
    }
}
write.csv(result_ratio, file = "num_selected_larger_than_1.csv", row.names = FALSE)

mean_num_select <- NULL
for (rw in c(120, 180)) {
    for (h in c(0.0834, 0.25, 0.5, 1, 2, 3)) {
        variable_name <- paste0("Result_", rw, "_", h, "_cv")
        assign("result_temp", eval(as.symbol(variable_name)))
        result_temp$Lasso$beta_hat
        mean_num_select <- rbind(mean_num_select, 
                              c(mean(rowSums(result_temp$ALasso$beta_hat != 0)),
                                mean(rowSums(result_temp$RepLasso$beta_hat != 0)),
                                mean(rowSums(result_temp$Lasso$beta_hat != 0)),
                                mean(rowSums(result_temp$Lasso_Std$beta_hat != 0))
                              ))
    }
}

write.csv(mean_num_select, file = "mean_num_select.csv", row.names = FALSE)


load("Result_obj_bic.RData")

result_ratio_bic <- NULL
for (rw in c(120, 180)) {
    for (h in c(0.0834, 0.25, 0.5, 1, 2, 3)) {
        variable_name <- paste0("Result_", rw, "_", h, "_bic")
        assign("result_temp", eval(as.symbol(variable_name)))
        result_ratio_bic <- rbind(result_ratio_bic, 
                              c(mean(rowSums(result_temp$ALasso$beta_hat != 0) > 1),
                                mean(rowSums(result_temp$RepLasso$beta_hat != 0) > 1),
                                mean(rowSums(result_temp$Lasso$beta_hat != 0) > 1),
                                mean(rowSums(result_temp$Lasso_Std$beta_hat != 0) > 1)
                              ))
    }
}
write.csv(result_ratio_bic, file = "num_selected_larger_than_1_bic.csv", row.names = FALSE)

mean_num_select_bic <- NULL
for (rw in c(120, 180)) {
    for (h in c(0.0834, 0.25, 0.5, 1, 2, 3)) {
        variable_name <- paste0("Result_", rw, "_", h, "_bic")
        assign("result_temp", eval(as.symbol(variable_name)))
        mean_num_select_bic <- rbind(mean_num_select_bic, 
                                 c(mean(rowSums(result_temp$ALasso$beta_hat != 0)),
                                   mean(rowSums(result_temp$RepLasso$beta_hat != 0)),
                                   mean(rowSums(result_temp$Lasso$beta_hat != 0)),
                                   mean(rowSums(result_temp$Lasso_Std$beta_hat != 0))
                                 ))
    }
}

write.csv(mean_num_select_bic, file = "mean_num_select_bic.csv", row.names = FALSE)
