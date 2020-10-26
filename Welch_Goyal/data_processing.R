# Data Processing
# Koo's data

library(readxl)
library(zoo)
library(lubridate)
library(tidyverse)


for(h in c(0.0834, 0.25, 0.5, 1:6)){
    m = as.integer(h*12)
    
    D = read.csv("Goyal_Koo.csv")
    D$yyyymm = as.yearmon( ymd(D$yyyymm, truncated = 2) )
    
    T = nrow(D)
    LongReturn = rep(NA, T)
    for(t in 1:(T-m+1)){
        LongReturn[t] =  sum( D$equity_premium[t:(t+m-1)] )
    }
    D$LongReturn = LongReturn
    
    Monthly.Data = D[-1, ]
    
    
    save(Monthly.Data, file = paste0("Monthly.koo.h", h, ".RData") )
    
}
