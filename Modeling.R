
########################################################################################
              # UCHICAGO - TIME-SERIES ANALYSIS
              # TOPIC: DROUGHT FORECASTING
              # Team Memeber: Nazih Kalo, Markus Wehr,Tam Nguyen, Vineela Potula
########################################################################################


########################################################################################
                                # IMPORT LIBRARY
########################################################################################

library(TSA)
library(fpp)
library(fpp2)
library(tseries)
library(zoo)
library(dplyr)
library(tidyverse)
library(gridExtra)
library(reshape2)
library(zoo)
library(xts)
library(ggplot2)
library(vars)
library(forecast)
library(tidyr)
library(knitr)
library(FinTS)
library(fGarch)

########################################################################################
                                # IMPORT DATA
########################################################################################

# Somalia SPEI index
Somalia<-read.csv('Somalia_Mogadishu_SPEI_long2.25_lat44.25.csv')
# Ethiopia SPEI index
Ethiopia <- read.csv('Ethiopia_AddisAbaba_SPEI_long9.75_lat38.75.csv')


# SPEI DATA
# Transform Somalia SPEI information
spei.Somalia<-Somalia[,1:2] #Restrict to the 1 month SPEI
spei.Ethiopia<-Ethiopia[,1:2] #Restrict to the 1 month SPEI


# Convert to TS object
spei.Somalia.ts <- ts(data = spei.Somalia[,2], start = c(1950,1), end = c(2018,11), frequency = 12)
spei.Ethiopia.ts <- ts(data = spei.Ethiopia[,2], start = c(1950,1), end = c(2018,11), frequency = 12)

# Somalia & Ethiopia SPEI Timeplot
autoplot(cbind(spei.Somalia.ts, spei.Ethiopia.ts), facets = TRUE, ylab = 'SPEI', main = 'Somalia & Ethiopia SPEI Timeplot') + 
  geom_vline(aes(xintercept = 2000, colour = 'Cutoff'), linetype = "dashed")


########################################################################################
                                  # Train Test Split
########################################################################################


# Train and test Somalia
train.Somalia <- window(spei.Somalia.ts, start = c(1999, 12), end = c(2016,11))
test.Somalia <-  window(spei.Somalia.ts, start =c(2016,12))


# Train and test Ethiopia
train.Ethiopia <- window(spei.Ethiopia.ts, start = c(1999, 12), end = c(2016,11))
test.Ethiopia <-  window(spei.Ethiopia.ts, start =c(2016,12))


# Show the train test split
tt_S <- autoplot(cbind(train.Somalia, test.Somalia), main = 'Somalia-Mogadishu Monthly SPEI', ylab = 'SPEI')
tt_E <- autoplot(cbind(train.Ethiopia, test.Ethiopia), main = 'Ethiopia-Addis Ababa Monthly SPEI', ylab = 'SPEI')
grid.arrange(tt_S, tt_E)


# No differencing required
ndiffs(train.Somalia); ndiffs(test.Somalia) 
nsdiffs(train.Somalia); nsdiffs(test.Somalia) 

ndiffs(train.Ethiopia); ndiffs(test.Ethiopia) 
nsdiffs(train.Ethiopia); nsdiffs(test.Ethiopia)

BoxCox.lambda(train.Somalia); BoxCox.lambda(test.Somalia) 
BoxCox.lambda(train.Ethiopia); BoxCox.lambda(test.Ethiopia) 

# ADF, KPSS test
adf.test(train.Ethiopia)
kpss.test(train.Ethiopia)

# Plot variance overtime
var.somalia <- autoplot(rollapply(cbind(train.Somalia,test.Somalia), width = 2, FUN = var), ylab = 'Variance', main = 'Somalia Variance Over Time')
var.ethiopia <- autoplot(rollapply(cbind(train.Ethiopia, test.Ethiopia), width = 2, FUN = var), ylab = 'Variance', main = 'Ethiopia Variance Over Time')
grid.arrange(var.somalia, var.ethiopia)

########################################################################################
                            # SEASONAL PLOTS 
#######################################################################################

ggseasonplot(train.Ethiopia, year.labels=TRUE, year.labels.left=TRUE) + ylab("SPEI") +
  ggtitle("Seasonal plot: Somalia SPEI_1 Value")

ggseasonplot(train.Ethiopia, polar=TRUE) +
  ylab("SPEI") +
  ggtitle("Polar seasonal plot: Somalia SPEI_1 Value")

ggsubseriesplot(train.Ethiopia) +
  ylab("SPEI") +
  ggtitle("Seasonal subseries plot: Somalia SPEI Value")


########################################################################################
                            # AUTOCORRELATION | ACF | PACF
########################################################################################

#Get ACF and PACF plot - SOMALIA
ggtsdisplay(train.Somalia, main = 'Somalia ACF/PACF', ylab = 'SPEI') 
ggAcf(spei.Somalia.ts, lag=96) #ACF plot for up to 96 lags (8 years)
ggPacf(spei.Somalia.ts, lag=96) #PACF plot for up to 96 lags (8 years)
#The ACF plot is decaying very slowly, this could be indicative of a ARIMA model. 


#Get ACF and PACF plot - ETHIOPIA
ggtsdisplay(train.Ethiopia,  main = 'Ethiopia ACF/PACF', ylab = 'SPEI') 
ggAcf(spei.Ethiopia.ts, lag=96) #ACF plot for up to 96 lags (8 years)
ggPacf(spei.Ethiopia.ts, lag=96) #PACF plot for up to 96 lags (8 years)


########################################################################################
                                  # Moving Averages
########################################################################################

autoplot(decompose(train.Somalia))
autoplot(cbind(ma(ma(train.Somalia, order = 12, centre = FALSE), order =2),
               ma(train.Somalia, order = 3), 
               train.Somalia), facets = TRUE)


###########################################################################################
                                # MODELING
                                # 1. BENCHMARK MODELS
###########################################################################################


###########################################
#SOMALIA

forecast_horizons <- c(1:24)
# SOMALIA FORECAST ACCURACY FOR ALL HORIZONS
benchmark_models_s <- list()
for (h in forecast_horizons) {
  naive_forc <- naive(train.Somalia, h = h)
  mean_forc <- meanf(train.Somalia, h = h)
  seas_naive_forc <- snaive(train.Somalia, h = h)
  naive_drift_forc <- rwf(train.Somalia, h = h, drift = TRUE)
  naive_acc <- accuracy(f = naive_forc, x = test.Somalia)[2,2]
  mean_acc <- accuracy(f = mean_forc, x = test.Somalia)[2,2]
  seas_naive_acc <- accuracy(f = seas_naive_forc, x = test.Somalia)[2,2]
  naive_drift_acc <- accuracy(f = naive_drift_forc, x = test.Somalia)[2,2]
  benchmark_models_s[[as.character(h)]] <- c(naive_acc,mean_acc,seas_naive_acc,naive_drift_acc)
}

meanf(train.Somalia, h = 24)
meanf(train.Ethiopia, h = 24)

# Combine all benchmark models together
benchmark_models_s.df <-  as.data.frame(benchmark_models_s, row.names = c('Naive', 'Mean', 'Seasonal Naive', 'Naive w/ Drift'))
colnames(benchmark_models_s.df) <- 1:24
bench.s.df = cbind(benchmark_models_s.df, Model = rownames(benchmark_models_s.df))
bench.s.df_melt = melt(bench.s.df, id = "Model")


# Plot RMSE Over Forecasting Horizons - Somalia
ggplot(bench.s.df_melt, aes(x = variable, y = value, group = Model)) + geom_line(aes(color=as.factor(Model))) + 
  ggtitle('RMSE Over Forecasting Horizons - Somalia')+
  xlab("Horizon") + ylab("RMSE") + scale_x_discrete('Horizon', seq(0,25,5) , seq(0,25,5))+ labs(color = "Model")


# Table RMSE for each model (12 & 24 months)
tbl <- tableGrob(round(benchmark_models_s.df[,c(12,24)],4), cols = c('h=12', 'h=24'), rows = c('Naive', 'Mean', 'Seasonal Naive', 'Naive w/ Drift'))
grid.arrange(tbl, as.table=TRUE)
tbl_means <- tableGrob(round(rowMeans(benchmark_models_s.df),4),cols = 'RMSE Overall' ,rows = c('Naive', 'Mean', 'Seasonal Naive', 'Naive w/ Drift'))
grid.arrange(tbl_means, as.table=TRUE)


# Plot Benchmark Forecasts for Somalia SPEI Value
bench_som <- autoplot(window(spei.Somalia.ts, start = c(2010, 12))) + autolayer(naive_forc, series = 'Naive', PI = FALSE) +
  autolayer(mean_forc, series = 'Mean',PI = FALSE) +
  autolayer(seas_naive_forc, series = 'Seasonal Naive',PI = FALSE) +
  autolayer(naive_drift_forc, series = 'Naive w/ Drift',PI = FALSE) +
  ggtitle("Benchmark Forecasts for Somalia SPEI Value") + xlab("Year") + ylab("SPEI") + guides(colour=guide_legend(title="Forecast"))


tbl <- tableGrob(round(benchmark_models_s.df,4), rows = c('Naive', 'Mean', 'Seasonal Naive', 'Naive w/ Drift'))
grid.arrange(bench_som, tbl,ncol = 1,
             nrow=2,
             as.table=TRUE, heights=c(3,1))


#####################################################
# ETHIOPIA

# ETHIOPIA FORECAST ACCURACY FOR ALL HORIZONS
benchmark_models_e <- list()
for (h in forecast_horizons) {
  naive_forc.e <- naive(train.Ethiopia, h = h)
  mean_forc.e <- meanf(train.Ethiopia, h = h)
  seas_naive_forc.e <- snaive(train.Ethiopia, h = h)
  naive_drift_forc.e <- rwf(train.Ethiopia, h = h, drift = TRUE)
  naive_acc.e <- accuracy(f = naive_forc.e, x = test.Ethiopia)[2,2]
  mean_acc.e <- accuracy(f = mean_forc.e, x = test.Ethiopia)[2,2]
  seas_naive_acc.e <- accuracy(f = seas_naive_forc.e, x = test.Ethiopia)[2,2]
  naive_drift_acc.e <- accuracy(f = naive_drift_forc.e, x = test.Ethiopia)[2,2]
  benchmark_models_e[[as.character(h)]] <- c(naive_acc.e,mean_acc.e,seas_naive_acc.e,naive_drift_acc.e)
}


# Combine all benchmark models together
benchmark_models_e.df <-  as.data.frame(benchmark_models_e, row.names = c('Naive', 'Mean', 'Seasonal Naive', 'Naive w/ Drift'))
colnames(benchmark_models_e.df) <- 1:24
bench.e.df = cbind(benchmark_models_e.df, Model = rownames(benchmark_models_e.df))
bench.e.df_melt = melt(bench.e.df, id = "Model")


# Plot RMSE Over Forecasting Horizons - Ethiopia
ggplot(bench.e.df_melt, aes(x = variable, y = value, group = Model)) + geom_line(aes(color=as.factor(Model))) + 
  ggtitle('RMSE Over Forecasting Horizons - Ethiopia')+
  xlab("Horizon") + ylab("RMSE") + scale_x_discrete('Horizon', seq(0,25,5) , seq(0,25,5))+ labs(color = "Model")


# Plot Benchmark Forecasts for Ethiopia SPEI Value
bench_eth <- autoplot(window(spei.Ethiopia.ts, start = c(2010, 12))) + autolayer(naive_forc.e, series = 'Naive', PI = FALSE) +
  autolayer(mean_forc.e, series = 'Mean',PI = FALSE) +
  autolayer(seas_naive_forc, series = 'Seasonal Naive',PI = FALSE) +
  autolayer(naive_drift_forc.e, series = 'Naive w/ Drift',PI = FALSE) +
  ggtitle("Benchmark Forecasts for Ethiopia SPEI Value") + xlab("Year") + ylab("SPEI") + guides(colour=guide_legend(title="Forecast"))


tbl <- tableGrob(round(benchmark_models_e.df[,c(12,24)],4), cols = c('h=12', 'h=24'), rows = c('Naive', 'Mean', 'Seasonal Naive', 'Naive w/ Drift'))
grid.arrange(tbl, as.table=TRUE)

tbl_means <- tableGrob(round(rowMeans(benchmark_models_e.df),4),cols = 'RMSE Overall' ,rows = c('Naive', 'Mean', 'Seasonal Naive', 'Naive w/ Drift'))
grid.arrange(tbl_means, as.table=TRUE)
grid.arrange(bench_som, bench_eth)



########################################################################################
                        # 2. Exponential Smoothing (SES, ETS)
# This method is suitable for forecasting data with no clear trend or seasonal pattern
########################################################################################


###########################
# ETHIOPIA

# SES_fit model
ses_E_12m <- ses(train.Ethiopia,h = 12)
ses_E_24m <- ses(train.Ethiopia,h = 24)
# SES_accuracy
accuracy(ses_E_12m,test.Ethiopia[1:12])
accuracy(ses_E_24m,test.Ethiopia[1:24])

# ETS_fit & forecast
ets_E_12m <- forecast(ets(train.Ethiopia,model=('ANN')), h=12)
ets_E_24m <- forecast(ets(train.Ethiopia,model=('ANN')), h=24)
# ETS_accuracy
accuracy(ets_E_12m,test.Ethiopia[1:12])
accuracy(ets_E_24m,test.Ethiopia[1:24])


# SES_Plot forecast Ethiopia
autoplot(ses_E_24m) +
  autolayer(fitted(ses_E_24m), series = "Fitted_24m") +
  ylab("SPEI_Ethiopia") + xlab("Year")

# ETS_Plot forecast Ethiopia
autoplot(ets_E_24m) + 
  autolayer(fitted(ets_E_24m), series = "Fitted_24m") +
  ylab("SPEI_Ethiopia") + xlab("Year")

# Plot - combine ETS and SES
autoplot(train.Ethiopia, xlab = 'Year',ylab = 'SPEI', main = 'Exponential Smoothing Forecast_Ethiopia') +
  autolayer(ses_E_24m$mean, series = "SES")+
  autolayer(ets_E_24m$mean, series = "ETS_ANN")


# SES_CV function 
model_ses_E <- function(x, h){ses(x,h = h)}
rmse_ses_ew_12m_E <- sqrt(colMeans((tsCV(train.Ethiopia, model_ses_E, h=12))^2,na.rm=TRUE))
rmse_ses_ew_24m_E <- sqrt(colMeans((tsCV(train.Ethiopia, model_ses_E, h=24))^2,na.rm=TRUE))
# ETS_function
model_ets_E <- function(x, h){forecast(ets(x,model=('ANN')), h=h)}
rmse_ets_ew_12m_E <- sqrt(colMeans((tsCV(train.Ethiopia, model_ets_E, h=12))^2,na.rm=TRUE))
rmse_ets_ew_24m_E <- sqrt(colMeans((tsCV(train.Ethiopia, model_ets_E, h=24))^2,na.rm=TRUE))

# RMSE table
E_ES_names <- c("ETS_ANN", "SES")
cv_rmse_df_E <- as.data.frame(cbind(rmse_ets_ew_24m_E,rmse_ses_ew_24m_E))
cv_rmse_df_E$Horizon <- 1:24
colnames(cv_rmse_df_E) <- c(E_ES_names, 'Horizon')

# plot cross validation
ggplot(data = melt(cv_rmse_df_E, id.vars = 'Horizon'), aes(x = Horizon, y = value, group = variable)) + 
  geom_line(aes(color = variable))+ xlab('Horizon') + ylab('CV RMSE') + ggtitle('CV RMSE over different forecast horizons_Ethiopia')+
  labs(color = "Model")

# plot Test RMSE over Horizon
ggplot(data = melt(cv_rmse_df_E_test, id.vars = 'Horizon'), aes(x = Horizon, y = value)) + 
  geom_line(aes())+ xlab('Test Horizon') + ylab('Test RMSE') + ggtitle('EST_ANN - Test RMSE over Horizon')

# Best Smooth Exponential Smoothing model for ETHIOPIA
RMSE.Ethiopia.best<-c()
forecast_horizon <- c(1:24)
for (h in forecast_horizon) {
  est_fc_E <- forecast(ets_E, h = h)
  ets_fc_E_acc <- accuracy(f = est_fc_E, x = test.Ethiopia,h = h)[2,2] 
  RMSE.Ethiopia.best<-c(RMSE.Ethiopia.best,ets_fc_E_acc)
}

# RMSE & AIcc table
Method<- c("SES", "ETS_ANN")
Ethiopia_ES <- data.frame(cbind(Method,rbind(
  round(rmse_ses_ew_24m_E[12],4),
  round(rmse_ets_ew_24m_E[12],4)),rbind(
    round(rmse_ses_ew_24m_E[24],4),
    round(rmse_ets_ew_24m_E[24],4)), rbind(
      round(ses_E_24m[["model"]][["aicc"]],4),
      round(ets_E_24m[["model"]][["aicc"]],4)
    )))
names(Ethiopia_ES)[2] <- "h = 12"
names(Ethiopia_ES)[3] <- "h = 24"
names(Ethiopia_ES)[4] <- "AICc"
Ethiopia_ES


####################
# SOMALIA #

# SES_fit model
ses_S_12m <- ses(train.Somalia,h = 12)
ses_S_24m <- ses(train.Somalia,h = 24)
# SES_accuracy
accuracy(ses_S_12m,test.Somalia[1:12])
accuracy(ses_S_24m,test.Somalia[1:24])

# ETS_fit & forecast
ets_S_12m <- forecast(ets(train.Somalia,model=('ANN')), h=12)
ets_S_24m <- forecast(ets(train.Somalia,model=('ANN')), h=24)
# ETS_accuracy
accuracy(ets_S_12m,test.Somalia[1:12])
accuracy(ets_S_24m,test.Somalia[1:24])


# SES_Plot forecast Somalia
autoplot(ses_S_24m) +
  autolayer(fitted(ses_S_24m), series = "Fitted_24m") +
  ylab("SPEI_Somalia") + xlab("Year") 

# ETS_Plot forecast Somalia
autoplot(ets_S_24m) + 
  autolayer(fitted(ets_S_24m), series = "Fitted_24m") +
  ylab("SPEI_Ethiopia") + xlab("Year")

# Plot - combine ETS and SES
autoplot(train.Somalia, xlab = 'Year',ylab = 'SPEI', main = 'Exponential Smoothing Forecast_Somalia') +
  autolayer(ses_S_24m$mean, series = "SES")+
  autolayer(ets_S_24m$mean, series = "ETS_ANN")


# CV function - SES
model_ses_S <- function(x, h){ses(x,h = h)}
rmse_ses_ew_12m_S <- sqrt(colMeans((tsCV(train.Somalia, model_ses_S, h=12))^2,na.rm=TRUE))
rmse_ses_ew_24m_S <- sqrt(colMeans((tsCV(train.Somalia, model_ses_S, h=24))^2,na.rm=TRUE))
# CV function - ETS
model_ets_S <- function(x, h){forecast(ets(x,model=('ANN')), h=h)}
rmse_ets_ew_12m_S <- sqrt(colMeans((tsCV(train.Somalia, model_ets_S, h=12))^2,na.rm=TRUE))
rmse_ets_ew_24m_S <- sqrt(colMeans((tsCV(train.Somalia, model_ets_S, h=24))^2,na.rm=TRUE))

# RMSE table
E_ES_names <- c("ETS_ANN", "SES")
cv_rmse_df_E <- as.data.frame(cbind(rmse_ets_ew_24m_E,rmse_ses_ew_24m_E))
cv_rmse_df_E$Horizon <- 1:24
colnames(cv_rmse_df_E) <- c(E_ES_names, 'Horizon')


# plot cross validation
ggplot(data = melt(cv_rmse_df_S, id.vars = 'Horizon'), aes(x = Horizon, y = value, group = variable)) + 
  geom_line(aes(color = variable))+ xlab('Horizon') + ylab('CV RMSE') + ggtitle('CV RMSE over different forecast horizons_Somalia')+
  labs(color = "Model")

# plot Test RMSE over Horizon
ggplot(data = melt(cv_rmse_df_E_test, id.vars = 'Horizon'), aes(x = Horizon, y = value)) + 
  geom_line(aes())+ xlab('Test Horizon') + ylab('Test RMSE') + ggtitle('EST_ANN - Test RMSE over Horizon')

# Best Smooth Exponential Smoothing model for SOMALIA
RMSE.Somalia.best<-c()
forecast_horizon <- c(1:24)
for (h in forecast_horizon) {
  est_fc_S <- forecast(ets_S, h = h)
  ets_fc_S_acc <- accuracy(f = est_fc_S, x = test.Somalia,h = h)[2,2] 
  RMSE.Somalia.best<-c(RMSE.Somalia.best,ets_fc_S_acc)
}

# RMSE & AIcc table
Method<- c("SES","ETS_ANN")
Somalia_ES <- data.frame(cbind(Method,rbind(
  round(rmse_ses_ew_24m_S[12],4),
  round(rmse_ets_ew_24m_S[12],4)),rbind(
    round(rmse_ses_ew_24m_S[24],4),
    round(rmse_ets_ew_24m_S[24],4)),rbind(
      round(ses_S_24m[["model"]][["aicc"]],4),
      round(ets_S_24m[["model"]][["aicc"]],4)
    )))
names(Somalia_ES)[2] <- "h = 12"
names(Somalia_ES)[3] <- "h = 24"
names(Somalia_ES)[4] <- "AICc"
Somalia_ES

###########################################################################################
                                      # 3. ARIMA, sARIMA
###########################################################################################

################
#3a. ARIMA
################

# SOMALIA #

# Automodels
auto.arrima_somalia <- auto.arima(train.Somalia, seasonal = TRUE, stepwise = FALSE, approximation = FALSE)
autoplot(auto.arrima_somalia, type = c('both'))

# Check ACF/PACF for model suggestions
ggtsdisplay(train.Somalia, main = 'Somalia Train Set', ylab = 'SPEI')
ggtsdisplay(train.Ethiopia, main = 'Ethiopia Train Set', ylab = 'SPEI')
typeof(eacf(train.Somalia)); eacf(train.Ethiopia)

ggAcf(train.Somalia, lag.max = 36,  main = 'Somalia Train Set ACF')
ggPacf(train.Somalia, lag.max = 36,  main = 'Somalia Train Set PACF')

ggAcf(train.Ethiopia, lag.max = 36,  main = 'Ethiopia Train Set ACF')
ggPacf(train.Ethiopia, lag.max = 36,  main = 'Ethiopia Train Set PACF')

as.data.frame(eacf(train.Somalia)$symbol)
tbl1 <- tableGrob(as.data.frame(eacf(train.Ethiopia)$symbol))
grid.arrange(ggPacf(train.Ethiopia, lag.max = 36), tbl1, 
             ncol =2, nrow =2, heights = c(5,1), widths = c(3,3),
             as.table=TRUE)


# check residuals
checkresiduals(train.Somalia)

# Manually Selected Models Somalia
ggAcf(train.Somalia)

# Try different order
ar.s.202 <- Arima(y = train.Somalia,order = c(2,0,2))
ar.s.103 <- Arima(y = train.Somalia,order = c(1,0,3))
ar.s.203 <- Arima(y = train.Somalia,order = c(2,0,3))
ar.s.302 <- Arima(y = train.Somalia,order = c(3,0,2))
ar.s.303 <- Arima(y = train.Somalia,order = c(3,0,3))
ar.s.102.100 <- Arima(y = train.Somalia,order = c(1,0,2), seasonal = list(order = c(1,0,1), period = 12))

#som_arima_models <- c('auto.arrima_somalia', 'ar.s.202','ar.s.103','ar.s.203','ar.s.302','ar.s.303', 'ar.s.102.100')
som_arima_names <- c('auto.ARIMA(1,0,2)', 'ARIMA(2,0,2)','ARIMA(1,0,3)','ARIMA(2,0,3)','ARIMA(3,0,2)','ARIMA(3,0,3)', 'ARIMA(1,0,2)(1,0,0)')
arima_som_aicc<- c(auto.arrima_somalia$aicc, ar.s.202$aicc, ar.s.103$aicc , ar.s.203$aicc , ar.s.302$aicc , ar.s.303$aicc, ar.s.102.100$aicc)
forecast_horizons <- c(1:24)

#SOMALIA FORECAST ACCURACY FOR ALL HORIZONS
somalia_models <- list()
for (h in forecast_horizons) {
  forc.1 <- forecast(auto.arrima_somalia, h = h)
  forc.2 <- forecast(ar.s.202, h = h)
  forc.3 <- forecast(ar.s.103, h = h)
  forc.4 <- forecast(ar.s.203, h = h)
  forc.5 <- forecast(ar.s.302, h = h)
  forc.6 <- forecast(ar.s.303, h = h)
  forc.7 <- forecast(ar.s.102.100, h = h)
  rmse.1 <- accuracy(f = forc.1, x = test.Somalia)[2,2]
  rmse.2 <- accuracy(f = forc.2, x = test.Somalia)[2,2]
  rmse.3 <- accuracy(f = forc.3, x = test.Somalia)[2,2]
  rmse.4 <- accuracy(f = forc.4, x = test.Somalia)[2,2]
  rmse.5 <- accuracy(f = forc.5, x = test.Somalia)[2,2]
  rmse.6 <- accuracy(f = forc.6, x = test.Somalia)[2,2]
  rmse.7 <- accuracy(f = forc.7, x = test.Somalia)[2,2]
  somalia_models[[as.character(h)]] <- c(rmse.1,rmse.2,rmse.3,rmse.4,rmse.5,rmse.6, rmse.7)
}
arima.rmse.s.df <-  as.data.frame(somalia_models, row.names = som_arima_names)



s.df = cbind(arima.rmse.s.df, model_name = rownames(arima.rmse.s.df))
s.df_melt = melt(s.df, id = "model_name")

# Plot RMSE Over Forecasting Horizons
ggplot(s.df_melt, aes(x = variable, y = value, group = as.factor(model_name) )) + 
  geom_line(aes(color=as.factor(model_name)), stat = "identity") + 
  ggtitle('RMSE Over Forecasting Horizons - Somalia')+
  xlab("Horizon") + ylab("RMSE")


# Plot all arima models together
arima_som_plot <- autoplot(window(train.Somalia, start = c(2010,1))) + autolayer(forc.1, series = 'auto.ARIMA(1,0,2)', PI = FALSE) +
  autolayer(forc.2, series = 'ARIMA(2,0,2)',PI = FALSE) +
  autolayer(forc.3, series = 'ARIMA(1,0,3)',PI = FALSE) +
  autolayer(forc.4, series = 'ARIMA(2,0,3)',PI = FALSE) +
  autolayer(forc.5, series = 'ARIMA(3,0,2)',PI = FALSE) +
  autolayer(forc.6, series = 'ARIMA(3,0,3)',PI = FALSE) +
  autolayer(forc.7, series = 'ARIMA(1,0,2)(1,0,0)',PI = FALSE) +
  ggtitle("ARIMA Forecasts for Somalia SPEI Value") + xlab("Year") + ylab("SPEI") + guides(colour=guide_legend(title="Model"))


s.df = cbind(arima.rmse.s.df[1:2], model_name = som_arima_names)
s.df_melt = melt(s.df, id = "model_name")

# Plot RMSE Over Forecasting Horizons - Somalia
RMSE_som_plot <- ggplot(s.df_melt, aes(x = variable, y = value, group = as.factor(model_name) )) + 
  geom_line(aes(color=as.factor(model_name)), stat = "identity") + 
  ggtitle('RMSE Over Forecasting Horizons - Somalia')+
  xlab("Horizon") + ylab("RMSE")+ theme(legend.position = "none")

som_arima_names <- c('auto.ARIMA(1,0,2)', 'ARIMA(2,0,2)','ARIMA(1,0,3)','ARIMA(2,0,3)','ARIMA(3,0,2)','ARIMA(3,0,3)', 'ARIMA(1,0,2)(1,0,0)')

tbl <- tableGrob(round(arima.rmse.s.df,4), rows = som_arima_names)
grid.arrange(tbl,
             as.table=TRUE)


#Analyzing selected model
auto.ARIMA.Somalia.102.f <-  function(x, h){forecast(Arima(x,order = c(1,0,2)),h =h)}
ARIMA.Somalia.202.f <- function(x, h){forecast(Arima(x,order = c(2,0,2)),h =h)}
ARIMA.Somalia.103.f <- function(x, h){forecast(Arima(x, order = c(1,0,3)),h =h)}
ARIMA.Somalia.203.f <- function(x, h){forecast(Arima(x, order = c(2,0,3)),h =h)}
ARIMA.Somalia.302.f <- function(x, h){forecast(Arima(x, order = c(3,0,2)),h =h)}
ARIMA.Somalia.303.f <- function(x, h){forecast(Arima(x, order = c(3,0,3)),h =h)}
ARIMA.Somalia.102.101.f <- function(x, h){forecast(Arima(x, order = c(1,0,2), seasonal =  list(order = c(1,0,1), period = 12)),h =h)}

# Cross validation
e.102 <-tsCV(train.Somalia, auto.ARIMA.Somalia.102.f, h = 24)
e.202<- tsCV(train.Somalia, ARIMA.Somalia.202.f, h = 24)
e.103<- tsCV(train.Somalia, ARIMA.Somalia.103.f, h = 24)
e.203<- tsCV(train.Somalia, ARIMA.Somalia.203.f, h = 24)
e.302<- tsCV(train.Somalia, ARIMA.Somalia.302.f, h = 24)
e.303 <- tsCV(train.Somalia, ARIMA.Somalia.303.f, h = 24)
e.102.101 <- tsCV(train.Somalia, ARIMA.Somalia.102.101.f, h = 24)

# RMSE
cbind(train.Somalia, test.Somalia)
rmse.102 <- sqrt(colMeans(e.102^2, na.rm = T))
rmse.202 <- sqrt(colMeans(e.202^2, na.rm = T))
rmse.103 <- sqrt(colMeans(e.103^2, na.rm = T))
rmse.203 <- sqrt(colMeans(e.203^2, na.rm = T))
rmse.302 <- sqrt(colMeans(e.302^2, na.rm = T))
rmse.303 <- sqrt(colMeans(e.303^2, na.rm = T))
rmse.102.101 <- sqrt(colMeans(e.102.101^2, na.rm = T))

# Corss validation
cv.rmse.df.s <- as.data.frame(cbind(rmse.102, rmse.202,rmse.103,rmse.203,rmse.302,rmse.303,rmse.102.101))
cv.rmse.df.s$Horizon <- 1:24
colnames(cv.rmse.df.s) <- c(som_arima_names, 'Horizon')

# Plot CV RMSE over different forecast horizons
ggplot(data = melt(cv.rmse.df.s, id.vars = 'Horizon'), aes(x = Horizon, y = value, group = variable)) + 
  geom_line(aes(color = variable))+ xlab('Horizon') + ylab('CV RMSE') + ggtitle('CV RMSE over different forecast horizons')+
  labs(color = "Model")


# Table of AICc
h12_h24_aicc_somalia <- cbind(t(cv.rmse.df.s)[1:7,c(12,24)], arima_som_aicc)
colnames(h12_h24_aicc_somalia) <- c('h = 12','h = 24', 'AICc')
tbl <- tableGrob(round(h12_h24_aicc_somalia,4))
grid.arrange(tbl,
             as.table=TRUE)


#BEST MODEL SOMALIA PLOT

autoplot(ar.s.303)

checkresiduals(ar.s.303)
autoplot(forecast(ar.s.303, h = 24))

best_somalia_arima <- as.data.frame(t(arima.rmse.s.df[6,]))

colnames(best_somalia_arima)
ggplot(best_somalia_arima, aes(x=1:24, y = best_somalia_arima$`ARIMA(3,0,3)`)) + geom_line() + ylab('Test RMSE') + xlab('Test Horizon')+
  ggtitle('Somalia ARIMA(3,0,3) Test RMSE over Horizon')

complete_som_ts <- window(spei.Somalia.ts, start = c(1999, 12))

som_arima_final <- Arima(complete_som_ts, order = c(3,0,3))

autoplot(forecast(som_arima_final, h = 24))

########################
# ETHIOPIA

auto.arrima_ethiopia <- auto.arima(train.Ethiopia, seasonal = TRUE, stepwise = FALSE, approximation = FALSE,trace = TRUE)
ggtsdisplay(train.Ethiopia, main = 'Ethiopia Train Set', ylab = 'SPEI')
eacf(train.Ethiopia)


#Manually Selected Models Ethiopia
ar.e.100.102 <- Arima(y = train.Ethiopia, order = c(1,0,0), seasonal = list(order = c(1,0,2), period = 9))
ar.e.001.102 <- Arima(y = train.Ethiopia, order = c(0,0,1), seasonal = list(order = c(1,0,2), period = 9))
ar.e.101.102 <- Arima(y = train.Ethiopia, order = c(1,0,1), seasonal = list(order = c(1,0,2), period = 9))
ar.e.002.102 <- Arima(y = train.Ethiopia, order = c(0,0,2), seasonal = list(order = c(1,0,2), period = 9))
ar.e.200.102 <- Arima(y = train.Ethiopia, order = c(2,0,0), seasonal = list(order = c(1,0,2), period = 9))
ar.e.000.100 <- Arima(y = train.Ethiopia, order = c(0,0,0), seasonal = list(order = c(1,0,0), period = 9))
eth_arima_names <- c('auto.ARIMA(2,0,3)', 'sARIMA(1,0,0)(1,0,2)','sARIMA(0,0,1)(1,0,2)','sARIMA(1,0,1)(1,0,2)',
                     'sARIMA(0,0,2)(1,0,2)','sARIMA(2,0,0)(1,0,2)', 'sARIMA(0,0,0)(1,0,0)')


arima_eth_aicc<- c(auto.arrima_ethiopia$aicc, ar.e.100.102$aicc ,ar.e.001.102$aicc,ar.e.101.102$aicc,
                   ar.e.002.102$aicc,ar.e.200.102$aicc,ar.e.000.100$aicc)


#Ethiopia FORECAST ACCURACY FOR ALL HORIZONS
ethiopia_models <- list()
for (h in forecast_horizons) {
  forc.1 <- forecast(auto.arrima_ethiopia, h = h)
  forc.2 <- forecast(ar.e.100.102, h = h)
  forc.3 <- forecast(ar.e.001.102, h = h)
  forc.4 <- forecast(ar.e.101.102, h = h)
  forc.5 <- forecast(ar.e.002.102, h = h)
  forc.6 <- forecast(ar.e.200.102, h = h)
  forc.7 <- forecast(ar.e.000.100, h = h)
  rmse.1 <- accuracy(f = forc.1, x = test.Ethiopia)[2,2]
  rmse.2 <- accuracy(f = forc.2, x = test.Ethiopia)[2,2]
  rmse.3 <- accuracy(f = forc.3, x = test.Ethiopia)[2,2]
  rmse.4 <- accuracy(f = forc.4, x = test.Ethiopia)[2,2]
  rmse.5 <- accuracy(f = forc.5, x = test.Ethiopia)[2,2]
  rmse.6 <- accuracy(f = forc.6, x = test.Ethiopia)[2,2]
  rmse.7 <- accuracy(f = forc.7, x = test.Ethiopia)[2,2]
  ethiopia_models[[as.character(h)]] <- c(rmse.1,rmse.2,rmse.3,rmse.4,rmse.5,rmse.6, rmse.7)
}
arima.rmse.e.df <-  as.data.frame(ethiopia_models, row.names = eth_arima_names)

autoplot(window(train.Ethiopia, start = c(2010,1))) + autolayer(forc.1, series = 'auto.ARIMA(2,0,3)', PI = FALSE) +
  autolayer(forc.2, series = 'sARIMA(1,0,0)(1,0,2)[9]',PI = FALSE) +
  autolayer(forc.3, series = 'sARIMA(0,0,1)(1,0,2)[9]',PI = FALSE) +
  autolayer(forc.4, series = 'sARIMA(1,0,1)(1,0,2)[9]',PI = FALSE) +
  autolayer(forc.5, series = 'sARIMA(0,0,2)(1,0,2)[9]',PI = FALSE) +
  autolayer(forc.6, series = 'sARIMA(2,0,0)(1,0,2)[9]',PI = FALSE) +
  autolayer(forc.7, series = 'sARIMA(0,0,0)(1,0,0)[9]',PI = FALSE) +
  ggtitle("ARIMA Forecasts for Ethiopia SPEI Value") + xlab("Year") + ylab("SPEI") + guides(colour=guide_legend(title="Model"))


e.df = cbind(arima.rmse.e.df[1:2], model_name = eth_arima_names)
e.df_melt = melt(e.df, id = "model_name")

ggplot(e.df_melt, aes(x = variable, y = value, group = as.factor(model_name) )) + 
  geom_line(aes(color=as.factor(model_name)), stat = "identity") + ggtitle('ARIMA RMSE Over Forecasting Horizons - Ethiopia')+
  xlab("Horizon") + ylab("RMSE") #+ theme(legend.position = "none")

eth_arima_names <- c('auto.ARIMA(2,0,3)', 'sARIMA(1,0,0)(1,0,2)[9]','sARIMA(0,0,1)(1,0,2)[9]','sARIMA(1,0,1)(1,0,2)[9]',
                     'sARIMA(0,0,2)(1,0,2)[9]','sARIMA(2,0,0)(1,0,2)[9]', 'sARIMA(0,0,0)(1,0,0)[9]')

tbl_e <- tableGrob(round(arima.rmse.e.df,4), rows = eth_arima_names)
grid.arrange(tbl_e,
             as.table=TRUE)


auto.ARIMA.Ethiopia.203.f <-  function(x, h){forecast(Arima(x,order = c(2,0,3)),h =h)}
ARIMA.Ethiopia.100.f <- function(x, h){forecast(Arima(x,order = c(1,0,0), seasonal = list(order = c(1,0,2), period = 9)),h =h)}
ARIMA.Ethiopia.001.f <- function(x, h){forecast(Arima(x,order = c(0,0,1), seasonal = list(order = c(1,0,2), period = 9)),h =h)}
ARIMA.Ethiopia.101.f <- function(x, h){forecast(Arima(x,order = c(1,0,1), seasonal = list(order = c(1,0,2), period = 9)),h =h)}
ARIMA.Ethiopia.002.f <- function(x, h){forecast(Arima(x,order = c(0,0,2), seasonal = list(order = c(1,0,2), period = 9)),h =h)}
ARIMA.Ethiopia.200.f <- function(x, h){forecast(Arima(x,order = c(2,0,0), seasonal = list(order = c(1,0,2), period = 9)),h =h)}
ARIMA.Ethiopia.000.f <- function(x, h){forecast(Arima(x,order = c(0,0,0), seasonal = list(order = c(1,0,0), period = 9)),h =h)}

e.203.e <-tsCV(train.Ethiopia, auto.ARIMA.Ethiopia.203.f, h = 24)
e.100.e <- tsCV(train.Ethiopia, ARIMA.Ethiopia.100.f, h = 24)
e.001.e <- tsCV(train.Ethiopia, ARIMA.Ethiopia.001.f, h = 24)
e.101.e <- tsCV(train.Ethiopia, ARIMA.Ethiopia.101.f, h = 24)
e.002.e <- tsCV(train.Ethiopia, ARIMA.Ethiopia.002.f, h = 24)
e.200.e <- tsCV(train.Ethiopia, ARIMA.Ethiopia.200.f, h = 24)
e.000.e <- tsCV(train.Ethiopia, ARIMA.Ethiopia.000.f, h = 24)


rmse.203.e <- sqrt(colMeans(e.203.e^2, na.rm = T))
rmse.100.e <- sqrt(colMeans(e.100.e^2, na.rm = T))
rmse.001.e <- sqrt(colMeans(e.001.e^2, na.rm = T))
rmse.101.e <- sqrt(colMeans(e.101.e^2, na.rm = T))
rmse.002.e <- sqrt(colMeans(e.002.e^2, na.rm = T))
rmse.200.e <- sqrt(colMeans(e.200.e^2, na.rm = T))
rmse.000.e <- sqrt(colMeans(e.000.e^2, na.rm = T))


cv.rmse.df.e <- as.data.frame(cbind(rmse.203.e, rmse.100.e,rmse.001.e,rmse.101.e,rmse.002.e,rmse.200.e,rmse.000.e))
cv.rmse.df.e$Horizon <- 1:24
colnames(cv.rmse.df.e) <- c(eth_arima_names, 'Horizon')

# CV RMSE over different forecast horizons
ggplot(data = melt(cv.rmse.df.e, id.vars = 'Horizon'), aes(x = Horizon, y = value, group = variable)) + 
  geom_line(aes(color = variable))+ xlab('Horizon') + ylab('CV RMSE') + ggtitle('CV RMSE over different forecast horizons')+
  labs(color = "Model")

h12_h24_aicc_ethiopia <- cbind(t(cv.rmse.df.e)[1:7,c(12,24)], arima_eth_aicc)
colnames(h12_h24_aicc_ethiopia) <- c('h = 12','h = 24', 'AICc')

tbl <- tableGrob(round(h12_h24_aicc_ethiopia,4))#, rows = som_arima_names)
grid.arrange(tbl,
             as.table=TRUE)

##########################################
# BEST MODEL ETHIOPIA PLOT

checkresiduals(ar.e.000.100)
autoplot(forecast(ar.e.000.100, h = 24))

autoplot(forecast(ar.e.000.100,h=24), series = 'auto.ARIMA(2,0,3)', PI = TRUE) + autolayer(window(spei.Ethiopia.ts, start = c(1999,12)))

best_eth_arima <- as.data.frame(t(arima.rmse.e.df[7,]))

# Plot Ethiopia sARIMA(0,0,0)(1,0,0) Test RMSE over Horizon
ggplot(best_eth_arima, aes(x=1:24, y = best_eth_arima$`sARIMA(0,0,0)(1,0,0)`)) + geom_line() + ylab('Test RMSE') + xlab('Test Horizon')+
  ggtitle('Ethiopia sARIMA(0,0,0)(1,0,0) Test RMSE over Horizon')

complete_eth_ts <- window(spei.Ethiopia.ts, start = c(1999, 12))

eth_arima_final <- Arima(complete_eth_ts, order = c(1,0,0), seasonal = list(order = c(1,0,2), period = 9))

autoplot(forecast(eth_arima_final, h = 24))


#####################################################################
                                # 4. ARIMAX
#####################################################################

arimax_som <- auto.arima(y = train.Somalia, xreg = window(somalia.fatalities.ts, end = c(2015,11)),approximation = FALSE, stepwise = FALSE)
autoplot(forecast(arimax_som, xreg = window(somalia.fatalities.ts, start = c(2015,12),  h = 36)))
accuracy(forecast(arimax_som, xreg = window(somalia.fatalities.ts, start = c(2015,12),  h = 36)), x = test.Somalia)


arimax_som <- auto.arima(y = train.Somalia, xreg = window(somalia.food.ts, end = c(2015,11)), approximation = FALSE, stepwise = FALSE)

autoplot(forecast(arimax_som, xreg = window(somalia.food.ts, start = c(2015,12),  h = 36)))
accuracy(forecast(arimax_som, xreg = window(somalia.food.ts, start = c(2015,12),  h = 36)), x = test.Somalia)

arimax_eth <- auto.arima(y = window(train.Ethiopia, start = c(2000,1)), xreg = window(ethiopia.food.ts, end = c(2015,11)), approximation = FALSE, stepwise = FALSE)

autoplot(forecast(arimax_eth, xreg = window(ethiopia.food.ts, start = c(2015,12),  h = 36)))
accuracy(forecast(arimax_eth, xreg = window(ethiopia.food.ts, start = c(2015,12),  h = 36)), x = test.Ethiopia)


############################################################################################
                      # 5. SPECTRAL ANALYSIS - Marmonic Regression
############################################################################################

# When there are long seasonal periods, a dynamic regression with Fourier terms is often better than other models we have considered in this book.
# harmonic regression approach where the seasonal pattern is modelled using Fourier terms with short-term time series dynamics handled by an ARMA error.


##############################################
# ETHIOPIA
plots <- list()
for (i in seq(6)) {
  fit <- auto.arima(train.Somalia, xreg = fourier(train.Somalia, K = i), seasonal = FALSE, lambda = 0)
  plots[[i]] <- autoplot(forecast(fit, xreg=fourier(train.Somalia, K=i, h=24))) + 
    xlab(paste("K=",i," AICC=",round(fit[["aicc"]],2))) +
    ylab("SPEI") + xlab('Time')
}
bestfit <- list(aicc = Inf)

for (i in seq(25)) {
  fit <- auto.arima(train.Somalia, xreg = fourier(train.Somalia, K = i), seasonal = FALSE, lambda = 0)
  if(fit$aicc < bestfit$aicc)
    bestfit <- fit
  else break;
}

summary(bestfit)

# autoplot
autoplot(bestfit$model)
plots

# check residuals
checkresiduals(fit)


##############################################
# SOMALIA 

plots <- list()
for (i in seq(6)) {
  fit <- auto.arima(train.Ethiopia, xreg = fourier(train.Ethiopia, K = i), seasonal = FALSE, lambda = 0)
  plots[[i]] <- autoplot(forecast(fit, xreg=fourier(train.Ethiopia, K=i, h=24))) + 
    xlab(paste("K=",i," AICC=",round(fit[["aicc"]],2))) +
    ylab("SPEI") + xlab('Time')
}

bestfit <- list(aicc = Inf)

for (i in seq(25)) {
  fit <- auto.arima(train.Ethiopia, xreg = fourier(train.Ethiopia, K = i), seasonal = FALSE, lambda = 0)
  if(fit$aicc < bestfit$aicc)
    bestfit <- fit
  else break;
}

summary(bestfit)

# autoplot
autoplot(bestfit$model)
plots

# check residuals
checkresiduals(fit)

########################################################################################
                               # 6. TBATS 
########################################################################################


####################################
# ETHIOPIA

tbats_eth <- tbats(train.Ethiopia)
accuracy(forecast(tbats_eth, h = 24), x = test.Ethiopia)
autoplot(forecast(tbats_eth, h = 24), main = 'Forecast from TBATS for Ethiopia', ylab = 'SPEI')
ethiopia_tbats <- matrix(nrow = 24, ncol = 1, dimnames = list(1:24, c('RMSE')))

# Combine forecast & accuracy
for (h in forecast_horizons) {
  forc.tbats <- forecast(tbats_eth, h = h)
  rmse.tbats <- accuracy(f = forc.tbats, x = test.Ethiopia)[2,2]
  ethiopia_tbats[h] <- c(rmse.tbats)
}

# Plot Test RMSE Over Forecasting Horizons - Ethiopia
ggplot(data = as.data.frame(ethiopia_tbats), aes(x = 1:24, y = RMSE)) + geom_line() + ggtitle('Test RMSE Over Forecasting Horizons - Ethiopia')+
  xlab('Horizon')
checkresiduals(tbats_eth) + ggtitle('Residuals from TBATS - ETHIOPIA')
plot(1:24, ethiopia_tbats, type = 'l')

####################################
# SOMALIA

tbats_som <- tbats(train.Somalia)
accuracy(forecast(tbats_som, h = 24), x = test.Somalia)
autoplot(forecast(tbats_som, h = 24),main = 'Forecast from TBATS for Somalia', ylab = 'SPEI')

# Combine forecast & accuracy
somlia_tbats <- matrix(nrow = 24, ncol = 1, dimnames = list(1:24, c('RMSE')))
for (h in forecast_horizons) {
  forc.tbats <- forecast(tbats_som, h = h)
  rmse.tbats <- accuracy(f = forc.tbats, x = test.Somalia)[2,2]
  somlia_tbats[h] <- c(rmse.tbats)
}


# Plot Test RMSE Over Forecasting Horizons - Somalia
ggplot(data = as.data.frame(somlia_tbats), aes(x = 1:24, y = RMSE)) + geom_line() + ggtitle('Test RMSE Over Forecasting Horizons - Somalia')+
  xlab('Horizon')

plot(1:24, somlia_tbats, type = 'l')
checkresiduals(tbats_som)


###########################################################################
                    #READ DATA FOR FINAL MODELS
###########################################################################

Somalia.temp<-read.csv('FINAL DATA/Somalia_Temperature.csv')
Somalia.temp<-Somalia.temp[,1:2]
Somalia.temp.ts<-ts(data = Somalia.temp[,2], start = c(1991,1), end = c(2018,12), frequency = 12)
autoplot(Somalia.temp.ts)

Ethiopia.temp<-read.csv('FINAL DATA/Ethiopia-weather-NEW.csv')
Ethiopia.temp<-Ethiopia.temp[,1:2]
Ethiopia.temp.ts<-ts(data = Ethiopia.temp[,2], start = c(1996,1), end = c(2019,8), frequency = 12)
autoplot(Ethiopia.temp.ts)

train.Somalia.temp<-window(Somalia.temp.ts, start = c(1999, 12), end = c(2016,11))
test.Somalia.temp<-window(Somalia.temp.ts, start = c(2016, 12), end = c(2018,11))

train.Ethiopia.temp<-window(Ethiopia.temp.ts, start = c(1999, 12), end = c(2016,11))
test.Ethiopia.temp<-window(Ethiopia.temp.ts, start = c(2016, 12), end = c(2018,11))

autoplot(cbind(train.Somalia.temp, test.Somalia.temp))
autoplot(cbind(train.Ethiopia.temp, test.Ethiopia.temp))


##############################################################################################
                                      # MODEL SELECTION
##############################################################################################

##################################
# SOMALIA

#Bench
mean_forc.s <- meanf(train.Somalia, h = 24)

#ARIMA
ar.s.303 <- Arima(y = train.Somalia,order = c(3,0,3))
arima_fc_24.s <- forecast(ar.s.303, h = 24)

#ETS
ets.somalia <- ets(train.Somalia,model=('ANN'))
ets_fc_24.s <- forecast(ets.somalia, h = 24)

#ARIMA errors
arima.lm.spei.temp.Somalia.auto<-auto.arima(train.Somalia, xreg = train.Somalia.temp, approximation = FALSE, stepwise = FALSE)
arima_err_fc_24.s <-forecast(arima.lm.spei.temp.Somalia.auto, xreg = snaive(train.Somalia.temp, h=24)$mean, h = 24)

#TBATS
tbats_som <- tbats(train.Somalia)
tbats_fc_24.s <- forecast(tbats_som, h =24)

#Spectral Analysis
spectral_k5_s <- auto.arima(train.Somalia,xreg=fourier(train.Somalia,3),seasonal = FALSE)
spectral_k5_fc_s <- forecast(spectral_k5_s, xreg = fourier(train.Somalia,3,24))

#VAR
variables_s <- cbind(spei = train.Somalia, temp = train.Somalia.temp)
VAR_som <- VAR(variables_s, p = 6, type = "both")
VAR_som_fc = forecast(VAR_som, h=24)
#rmse_s =  accuracy(fc_b_s [["forecast"]][["spei"]][["mean"]], x= test_ethiopia)[1,2]

# Plot Best Model Forecasts for Somalia SPEI Value
autoplot(window(spei.Somalia.ts, start = c(2013, 12))) + autolayer(mean_forc.s, series = 'Mean', PI = FALSE) +
  autolayer(arima_fc_24.s, series = 'ARIMA',PI = FALSE) +
  autolayer(ets_fc_24.s, series = 'ETS',PI = FALSE) +
  autolayer(arima_err_fc_24.s, series = 'ARIMA errors',PI = FALSE) +
  autolayer(tbats_fc_24.s, series = 'TBATS',PI = FALSE) +
  autolayer(spectral_k5_fc_s, series = 'Spectral',PI = FALSE)+
  autolayer(VAR_som_fc$forecast$spei, series = 'VAR',PI = FALSE)+
  ggtitle("Best Model Forecasts for Somalia SPEI Value") + xlab("Year") + ylab("SPEI") + guides(colour=guide_legend(title="Forecast"))

modelnames_s <- c('Mean', 'ARIMA', 'ETS', 'ARIMA errors', 'TBATS', 'Spectral','VAR')
final_models_somalia <- matrix(nrow = 24, ncol = 7, dimnames = list(c(1:24), modelnames_s))
forecast_horizons <- 1:24

# forecast_horizons
for (h in forecast_horizons) {
  mean_forc <- meanf(train.Somalia, h = h)
  arima_fc_s <- forecast(ar.s.303, h = h)
  ets_fc_s <- forecast(ets.somalia, h = h)
  arima_err_fc_s <- forecast(arima.lm.spei.temp.Somalia.auto, xreg = snaive(train.Somalia.temp, h=h)$mean, h = h)
  tbats_fc_s <- forecast(tbats_som, h = h)
  spectral_fc_s <- forecast(spectral_k5_s,xreg = fourier(train.Somalia,3,h) ,h = h)
  VAR_fc_s <- forecast(VAR_som, h = h)$forecast$spei
  mean_acc <- accuracy(f = mean_forc, x = test.Somalia)[2,2]
  arima_fc_s_acc <- accuracy(f = arima_fc_s, x = test.Somalia)[2,2]
  ets_fc_s_acc <- accuracy(f = ets_fc_s, x = test.Somalia)[2,2]
  arima_err_fc_s_acc<- accuracy(f = arima_err_fc_s, x = test.Somalia)[2,2]
  tbats_fc_s_acc <- accuracy(f = tbats_fc_s, x = test.Somalia)[2,2]
  spectral_fc_s_acc <- accuracy(f = spectral_fc_s, x = test.Somalia)[2,2]
  VAR_fc_s_acc <- accuracy(f = VAR_fc_s, x = test.Somalia)[2,2]
  final_models_somalia[h,1] <- mean_acc
  final_models_somalia[h,2] <- arima_fc_s_acc
  final_models_somalia[h,3] <- ets_fc_s_acc
  final_models_somalia[h,4] <- arima_err_fc_s_acc
  final_models_somalia[h,5] <- tbats_fc_s_acc
  final_models_somalia[h,6] <- spectral_fc_s_acc
  final_models_somalia[h,7] <- VAR_fc_s_acc
  
}

rmse_12_24_s <- t(final_models_somalia[c(12,24),])
cbind(rmse_12_24_s, 'Mean' = colMeans(final_models_somalia))

tbl <- tableGrob(round(cbind(rmse_12_24_s, 'Mean' = colMeans(final_models_somalia)),4), cols = c('h=12', 'h=24', 'Mean'))
grid.arrange(tbl, as.table=TRUE)
final_models_somalia_melt <- melt(final_models_somalia, varnames = c('Horizon', 'Model'))

# Plot RMSE Over Forecasting Horizons - Somalia
ggplot(final_models_somalia_melt, aes(x = Horizon, y = value, group = Model)) + geom_line(aes(color=as.factor(Model))) + 
  ggtitle('RMSE Over Forecasting Horizons - Somalia')+
  xlab("Horizon") + ylab("RMSE")+ labs(color = "Model")

########################################
# FINAL PREDICTION SOMALIA

arima_final_som <- Arima(y = window(spei.Somalia.ts, start = c(1999,12)),order = c(3,0,3))
arima_final_som_fc <- forecast(arima_final_som, h = 24)

summer_pred_som <- as.data.frame(arima_final_som_fc)[c(9,21),]
tbl_summer <- tableGrob(round(summer_pred_som,4))#, cols = c('h=12', 'h=24', 'Mean'))
grid.arrange(tbl_summer, as.table=TRUE)

# Plot Final ARIMA(3,0,3) Prediction for 2019-2020
autoplot(arima_final_som_fc, main = 'Final ARIMA(3,0,3) Prediction for 2019-2020') + ylab('SPEI')

# Check residuals
checkresiduals(arima_final_som)

####################################################
# ETHIOPIA

#Bench
naive_forc.e <- naive(train.Ethiopia, h = 24)

#ARIMA
ar.e.000.100 <- Arima(y = train.Ethiopia, order = c(0,0,0), seasonal = list(order = c(1,0,0), period = 9))
arima_fc_24.e <- forecast(ar.e.000.100, h = 24)

#ETS
ets.ethiopia<- ets(train.Ethiopia,model=('ANN'))
ets_fc_24.e <- forecast(ets.ethiopia, h = 24)

#ARIMA ERRORS
arima.lm.spei.temp.Ethiopia.202<-Arima(train.Ethiopia, order = c(0,0,0), list(order = c(1,0,0), period = 9), xreg = train.Ethiopia.temp)
arima_err_fc_24.e <-forecast(arima.lm.spei.temp.Ethiopia.202, xreg = snaive(train.Ethiopia.temp, h=24)$mean, h = 24)

#TBATS
tbats_eth <- tbats(train.Ethiopia)
tbats_fc_24.e <- forecast(tbats_eth, h =24)

#Spectral Analysis
spectral_k4_e <- auto.arima(train.Ethiopia,xreg=fourier(train.Ethiopia,2),seasonal = FALSE)
spectral_k4_fc_e <- forecast(spectral_k4_e, xreg = fourier(train.Ethiopia,2,24))

#VAR
variables_e <- cbind(spei = train.Ethiopia, temp = train.Ethiopia.temp)
VAR_eth <- VAR(variables_e, p = 8, type = "both")
VAR_eth_fc = forecast(VAR_eth, h=24)

# Plot Best Model Forecasts for Ethiopia SPEI Value
autoplot(window(spei.Ethiopia.ts, start = c(2013, 12))) + autolayer(naive_forc.e, series = 'Naive', PI = FALSE) +
  autolayer(arima_fc_24.e, series = 'ARIMA',PI = FALSE) +
  autolayer(ets_fc_24.e, series = 'ETS',PI = FALSE) +
  autolayer(arima_err_fc_24.e, series = 'ARIMA errors',PI = FALSE) +
  autolayer(tbats_fc_24.e, series = 'TBATS',PI = FALSE) +
  autolayer(spectral_k4_fc_e, series = 'Spectral',PI = FALSE) +
  autolayer(VAR_eth_fc$forecast$spei, series = 'VAR',PI = FALSE) +
  ggtitle("Best Model Forecasts for Ethiopia SPEI Value") + xlab("Year") + ylab("SPEI") + guides(colour=guide_legend(title="Forecast"))


modelnames_e <- c('Naive', 'ARIMA', 'ETS', 'ARIMA errors', 'TBATS','Spectral','VAR')
final_models_ethiopia <- matrix(nrow = 24, ncol = 7, dimnames = list(c(1:24), modelnames_e))
forecast_horizons <- 1:24

# forecast_horizons
for (h in forecast_horizons) {
  naive_fc.e <-  naive(train.Ethiopia, h = h)
  arima_fc_e <- forecast(ar.e.000.100, h = h)
  ets_fc_e <- forecast(ets.ethiopia, h = h)
  arima_err_fc_e <- forecast(arima.lm.spei.temp.Ethiopia.202, xreg = snaive(train.Ethiopia.temp, h=h)$mean, h = h)
  tbats_fc_e <- forecast(tbats_eth, h = h)
  spectral_fc_e <- forecast(spectral_k4_e,xreg = fourier(train.Ethiopia,2,h) ,h = h)
  VAR_fc_e <- forecast(VAR_eth, h = h)$forecast$spei
  naive_fc_acc_e <- accuracy(f = naive_fc.e, x = test.Ethiopia)[2,2]
  arima_fc_e_acc <- accuracy(f = arima_fc_e, x = test.Ethiopia)[2,2]
  ets_fc_e_acc <- accuracy(f = ets_fc_e, x = test.Ethiopia)[2,2]
  arima_err_fc_e_acc<- accuracy(f = arima_err_fc_e, x = test.Ethiopia)[2,2]
  tbats_fc_e_acc <- accuracy(f = tbats_fc_e, x = test.Ethiopia)[2,2]
  spectral_fc_e_acc <- accuracy(f = spectral_fc_e, x = test.Ethiopia)[2,2]
  VAR_fc_e_acc <- accuracy(f = VAR_fc_e, x = test.Ethiopia)[2,2]
  final_models_ethiopia[h,1] <- naive_fc_acc_e
  final_models_ethiopia[h,2] <- arima_fc_e_acc
  final_models_ethiopia[h,3] <- ets_fc_e_acc
  final_models_ethiopia[h,4] <- arima_err_fc_e_acc
  final_models_ethiopia[h,5] <- tbats_fc_e_acc
  final_models_ethiopia[h,6] <- spectral_fc_e_acc
  final_models_ethiopia[h,7] <- VAR_fc_e_acc
}

rmse_12_24_e <- t(final_models_ethiopia[c(12,24),])

cbind(rmse_12_24_e, 'Mean' = colMeans(final_models_ethiopia))
tbl <- tableGrob(round(cbind(rmse_12_24_e, 'Mean' = colMeans(final_models_ethiopia)),4), cols = c('h=12', 'h=24', 'Mean'))#, rows = c('Naive', 'Mean', 'Seasonal Naive', 'Naive w/ Drift'))
grid.arrange(tbl, as.table=TRUE)
final_models_ethiopia_melt <- melt(final_models_ethiopia, varnames = c('Horizon', 'Model'))

# Plot RMSE Over Forecasting Horizons - Ethiopia
ggplot(final_models_ethiopia_melt, aes(x = Horizon, y = value, group = Model)) + geom_line(aes(color=as.factor(Model))) + 
  ggtitle('RMSE Over Forecasting Horizons - Ethiopia')+
  xlab("Horizon") + ylab("RMSE")+ labs(color = "Model")


########################################
# FINAL PREDICTION ETHIOPIA

tbats_eth_final <- tbats(window(spei.Ethiopia.ts, start = c(1999,12)))
tbats_eth_final_fc <- forecast(tbats_eth_final, h =24)

# Plot Final TBATS Prediction for 2019-2020
autoplot(tbats_eth_final_fc, main = 'Final TBATS Prediction for 2019-2020') + ylab('SPEI')
checkresiduals(tbats_eth_final)

summer_pred_eth <- as.data.frame(tbats_eth_final_fc)[c(9,21),]
tbl_summer <- tableGrob(round(summer_pred_eth,4))#, cols = c('h=12', 'h=24', 'Mean'))
grid.arrange(tbl_summer, as.table=TRUE)


####################################################################################
                              # ARCH/ GARCH
####################################################################################


###############################
# ETHIOPIA

## Tranform data to positive series
test_E_tf <- test.Ethiopia +4
train_E_tf <- train.Ethiopia +4

# Compute thr return of Ethiopia
ln_E<-log(train_E_tf)
return_E <-diff(ln_E)

# Ljung Box test shows the return is not white noise
Box.test(return_E, lag = 1, type = "Ljung")

# acf test
acf(return_E, main = "Return - Ethiopia")
tsdisplay(return_E)

# Plot return of SPEI
autoplot(return_E, ylab = 'Return - r', main = 'Plot of return of SPEI _ Ethiopia')

# fit arima model
fit1<-arima(return_E, order = c(1,0,0))
summary(fit1)
res1<-resid(fit1)
# variance prediction
predict(arima(return_E, order = c(1,0,0)),n.ahead = 24)

# fit garch model
model1 <- garch(x = return_E, order=c(1,1))
summary(model1)

# predict 24months ahead
arch.fit_E <- garchFit(~garch(1,1), data = return_E, trace = F)
summary(arch.fit_E)
garch_model_E <- predict(arch.fit_E,n.ahead = 24)
garch_model_E 

garch_model_E$Variance <- garch_model_E$standardDeviation^2
garch_model_E$RMSE <- sqrt((garch_model_E$meanError)^2)

# Variance in train set
var(return_E)
# Variance in test set
var(return_E_test)

# Variance in testing set - ETHIOPIA
Dataset_Ethiopia <- c("Training", "Test")
Ehiopia_Variance <- data.frame(cbind(Dataset_Ethiopia,(rbind(round(var(return_E),4),round(var(return_E_test),4)))))
names(Ehiopia_Variance)[2] <- "Variance"
Ehiopia_Variance


####################################
# SOMALIA

## Tranform data to positive series
test_S_tf <- test.Somalia +4
train_S_tf <- train.Somalia +4

# Compute thr return of Somalia
ln_S_test <-log(test_S_tf)
return_S_test <-diff(ln_S_test)

# Ljung Box test shows the return is not white noise
Box.test(return_S, lag = 1, type = "Ljung")

# acf test
acf(return_S, main = "Return - Somalia")
tsdisplay(return_S)

# Plot return of SPEI
autoplot(return_S, ylab = 'Return - r', main = 'Plot of return of SPEI _ Somalia')

# fit arima model
fit2<-arima(return_S, order = c(1,0,0), method ="CSS")
summary(fit2)
res2<-resid(fit2)
# variance prediction
predict(arima(return_S, order = c(1,0,0)),n.ahead = 24)

# fit garch model
model2 <- garch( x = return_S, order=c(1,1))
summary(model2)

# predict 24months ahead
arch.fit_S <- garchFit(~garch(1,1), data = return_S, trace = F)
summary(arch.fit_S)
garch_model_S <- predict(arch.fit_S,n.ahead = 24)
garch_model_S 

garch_model_S$Variance <- garch_model_S$standardDeviation^2
garch_model_S$RMSE <- sqrt((garch_model_S$meanError)^2)

# Variance in train set
var(return_S)
# Variance in test set
var(return_S_test)

# Variance in testing set - SOMALIA
Dataset_Somalia <- c("Training", "Test")
Ehiopia_Variance <- data.frame(cbind(Dataset_Somalia,(rbind(round(var(return_S),4),round(var(return_S_test),4)))))
names(Somalia_Variance)[2] <- "Variance"
Somalia_Variance

###########################################################




