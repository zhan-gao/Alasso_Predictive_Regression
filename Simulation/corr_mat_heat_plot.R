library(ggcorrplot)
library(ggplot2)
library(corrplot)
library(RColorBrewer)

# -----------------------
cov_mat_1 <- dgp_1_parameter$sigma_mat
colnames(cov_mat_1) <- rownames(cov_mat_1) <- c(paste0("e", 1:9), "u")

cov_mat_2 <- dgp_2_parameter$sigma_mat
colnames(cov_mat_2) <- rownames(cov_mat_2) <- c(paste0("z", 1:3), paste0("v", 1:4), paste0("e", 1:5) ,"u")
# -----------------------


heat_map_dgp_1 <- corrplot(cov2cor(cov_mat_1), 
                           method = "color")
heat_map_dgp_2 <- corrplot(cov2cor(cov_mat_2), 
                           method = "color")

ggsave(filename = "heat_map_dgp1.png",
       plot = heat_map_dgp_1,
       device = "png",
       width = 150,
       height = 150,
       unit = "mm")

ggsave(filename = "heat_map_dgp2.png",
       plot = heat_map_dgp_2,
       device = "png",
       width = 150,
       height = 150,
       unit = "mm")
