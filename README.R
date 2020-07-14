# ************************************************************************************
# ************************************************************************************
# Author: Gianluca Broll, University of Trento
# Date: 2020-07-09
# ************************************************************************************
# HOW DID COVID-19 CRISIS AFFECT FINANCIAL RISK MODELING?
# ************************************************************************************
# ************************************************************************************

remove(list=ls())
library(quantmod)
library(xts)
library(PerformanceAnalytics)
library(tseries)
library(rugarch)
library(ismev)
library(evir)
library(ggplot2)
library(grid)
library(gridExtra)

# ************************************************************************************
# RETRIEVING DATA
# ************************************************************************************

getSymbols("^GSPC", from = "1990-01-01", to = "2020-07-09", src =  "yahoo")
SP500            <- Return.calculate(GSPC$GSPC.Adjusted, method = "log")[-1,] * 100
SP500_preCov     <- SP500["/201912"]
SP500_postCov    <- SP500["199007/"]

ggplot(data = SP500/100, aes(x = index(SP500))) +
  geom_line(aes(y = GSPC.Adjusted), color = "lightslategray") +
  ggtitle("S&P 500", subtitle = "Daily returns from 1990-01-01 to 2020-07-09") +
  theme_classic() +
  theme(plot.title = element_text(size = 18, face = "bold")) +
  labs(caption = "Data source: Yahoo! Finance") +
  theme(plot.caption = element_text(face = "italic")) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Returns") + 
  xlab("Date")

garch             = ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(2,2)))
mdl               <- ugarchfit(spec = garch, data = SP500/100, solver = "hybrid")
std_res           <- residuals(mdl)/sigma(mdl)
data              <- data.frame(index(std_res), std_res)
colnames(data)    <- c("Date", "Residuals")

ggplot(data, aes(x = Date)) +
  geom_line(aes(y = Residuals), color = "lightslategray") +
  ggtitle("Standardized Residuals", subtitle = "SP500 from 1990-01-01 to 2020-07-09") +
  theme_classic() +
  theme(plot.title = element_text(size = 18, face = "bold")) +
  labs(caption = "Data source: Yahoo! Finance") +
  theme(plot.caption = element_text(face = "italic"))

# ************************************************************************************
# TAIL MODELLING (GPD ESTIMATION) 
# ************************************************************************************

threshold        <- 3
alpha            <- 0.999

# ************************************************************************************
# Pre Covid-19 crisis
# ************************************************************************************

garch_preCov     = ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(2,2)))
mdl_preCov       <- ugarchfit(spec = garch_preCov, data = SP500_preCov, solver = "hybrid")
residual_preCov  <- residuals(mdl_preCov)
cond_vol_preCov  <- sigma(mdl_preCov)
std_res_preCov   <- residual_preCov/cond_vol_preCov
data             <- data.frame(index(residual_preCov), std_res_preCov)
colnames(data)   <- c("Date", "Residuals")
value_preCov     <- as.vector(-std_res_preCov)
g_preCov         <- gpd.fit(value_preCov, threshold)
h                <- hist(g_preCov$data, plot = F)
x                <- seq(threshold, max(h$breaks), length = 100)
y                <- dgpd(x, xi = g_preCov$mle[2], beta = g_preCov$mle[1], mu = threshold)
q_Z_preCov       <- threshold + (g_preCov$mle[1] / g_preCov$mle[2]) * (((1 - alpha) / g_preCov$rate) ^ (-g_preCov$mle[2]) - 1)
data             <- data.frame(g_preCov$data)
colnames(data)   <- "Exceedances"
# xx               <- seq(threshold,8, by = 0.001)
# tt               <- truncnorm::dtruncnorm(xx, a = threshold, mean = mean(SP500_preCov), sd = sd(SP500_preCov))  # Uncomment to see the diffrence with normal

p1 <- ggplot(data) +
  geom_histogram(aes(x = Exceedances, y = ..density..), bins = 14, color = "lightslategray", fill = "lightslategray", alpha = 0.8) +
  geom_line(data = data.frame(x, y), aes(x, y), color = "darkblue") +
  xlab("Standardized Residual Exceedances") + 
  ylab("Density") +
  ggtitle("Fitted GPD to Standardized Residuals", subtitle = paste("Pre Covid-19 Crisis ", "(Shape = ", round(g_preCov$mle[2], 4), ", Scale = ", round(g_preCov$mle[1],4), ")", sep = "")) +
  theme_classic() +
  theme(plot.title = element_text(size = 18, face = "bold")) +
  geom_segment(aes(x = q_Z_preCov , y = 0, xend = q_Z_preCov, yend = 1, color = "q_Z_preCov"), linetype = "dashed") +
  scale_color_manual(name = "", labels = paste(alpha*100, "% Quantile", sep = ""), values = "red") +
  geom_label(aes(x = q_Z_preCov, y = 1.05, label = round(q_Z_preCov, 4)), color = "red") +
  theme(legend.position = "bottom")
# + geom_line(data = data.frame(xx, tt), aes(x = xx, y = tt))  # Uncomment to see the diffrence with normal

# ************************************************************************************
# Post Covid-19 Crisis
# ************************************************************************************

garch_postCov      = ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(2,2)))
mdl_postCov        <- ugarchfit(spec = garch_postCov, data = SP500_postCov, solver = "hybrid")
residual_postCov   <- residuals(mdl_postCov)
cond_vol_postCov   <- sigma(mdl_postCov)
std_res_postCov    <- residual_postCov/cond_vol_postCov
data               <- data.frame(index(residual_postCov), std_res_postCov)
colnames(data)     <- c("Date", "Residuals")
value_postCov      <- as.vector(-std_res_postCov)
g_postCov          <- gpd.fit(value_postCov, threshold)
h                  <- hist(g_postCov$data, plot = F)
x                  <- seq(threshold, max(h$breaks), length = 100)
y                  <- dgpd(x, xi = g_postCov$mle[2], beta = g_postCov$mle[1], mu = threshold)
q_Z_postCov        <- threshold + g_postCov$mle[1] / g_postCov$mle[2] * (((1 - alpha) / g_postCov$rate) ^ (-g_postCov$mle[2]) - 1)
data               <- data.frame(g_postCov$data)
colnames(data)     <- "Exceedances"
# xx                 <- seq(threshold,8, by = 0.001)
# tt                 <- truncnorm::dtruncnorm(xx, a = threshold, mean = mean(SP500_postCov), sd = sd(SP500_postCov))  # Uncomment to see the diffrence with normal

p2 <- ggplot(data) +
  geom_histogram(aes(x = Exceedances, y = ..density..), bins = 14, col = "lightslategray", fill = "lightslategray", alpha = 0.8) +
  geom_line(data = data.frame(x, y), aes(x, y), color  = "darkblue") +
  xlab("Standardized Residual Exceedances") + 
  ylab("Density") +
  theme_classic() +
  ggtitle("Fitted GPD to Standardized Residuals", subtitle = (paste("Post Covid-19 Crisis ", "(Shape = ", round(g_postCov$mle[2], 4), ", Scale = ", round(g_postCov$mle[1],4), ")", sep = ""))) +
  geom_segment(aes(x = q_Z_postCov , y = 0, xend = q_Z_postCov, yend = 1, color = "q_Z_postCov"), linetype = "dashed") +
  scale_color_manual(name = "", labels = paste(alpha*100, "% Quantile", sep = ""), values = "red") +
  theme(plot.title = element_text(size = 18, face = "bold")) +
  geom_label(aes(x = q_Z_postCov, y = 1.05, label = round(q_Z_postCov, 4)), color = "red") +
  theme(legend.position = "bottom")
  # + geom_line(data = data.frame(xx, tt), aes(x = xx, y = tt))      # Uncomment to see the diffrence with normal
  
grid.arrange(p1, p2, nrow = 1)

# ************************************************************************************
# RISK MEASURE ESTIMATE (VALUE AT RISK) 
# ************************************************************************************

threshold       <- 2
alpha           <- 0.99
start_date      <- as.numeric(which(index(SP500) == "2018-12-31"))

# ************************************************************************************
# VaR Forecast on a Rolling Window (250 obs), Forecasting Starting Date: 2018-12-31
# ************************************************************************************

garch_var       = ugarchspec(variance.model = list(model = "eGARCH"))
mdl_var         <- ugarchfit(spec = garch_var, data = SP500["/2018-12-31"]/100, solver = "hybrid")
std_res_var     <- residuals(mdl_var)/sigma(mdl_var)
g_var           <- gpd.fit(-std_res_var, threshold)
q_Z_var         <- threshold + g_var$mle[1] / g_var$mle[2] * (((1 - alpha) / g_var$rate) ^ (-g_var$mle[2]) - 1)
window          <- 250
VaR             <- as.vector(length(SP500["20181231/"]))

for (i in 1:(length(SP500["20181231/"])-1)) {
  data        <- SP500[(start_date - window + i - 1):(start_date + i - 1)]/100
  garch       <- ugarchspec(variance.model = list(model = "eGARCH"))
  mdl         <- ugarchfit(garch, data, solver = "hybrid")
  forecast    <- ugarchforecast(mdl, n.ahead = 1)
  
  if (i %% 30 == 0) {
    data_gpd = SP500[1:(start_date + i)]
    garch_gpd     <- ugarchspec(variance.model = list(model = "eGARCH"))
    mdl_gpd       <- ugarchfit(garch_gpd, data_gpd, solver = "hybrid")
    gpd           <- gpd.fit(-residuals(mdl_gpd)/sigma(mdl_gpd), threshold)
    q_Z_var       <- threshold + gpd$mle[1] / gpd$mle[2] * (((1 - alpha) / gpd$rate) ^ (-gpd$mle[2]) - 1)
  }
  
  VaR[i] <- as.numeric(forecast@forecast$seriesFor) +  as.numeric(forecast@forecast$sigmaFor) * q_Z_var
}

VaR_forecast           <- data.frame(index(SP500["20190102/"]), VaR, data = SP500["20190102/"])
colnames(VaR_forecast) <- c("Date", "VaR", "SP500")

ggplot(VaR_forecast, aes(x = Date)) +
  geom_col(aes(y = SP500/100), width = 1, color = "lightslategray", fill = "lightslategray") +
  geom_line(aes(y = -VaR), color = "darkblue") +
  scale_y_continuous(labels = scales::percent) +
  xlab("") +
  theme_classic() +
  ggtitle("Daily VaR Forecast S&P 500", subtitle = "Rolling Window on 250 Observations") +
  theme(plot.title = element_text(size = 18, face = "bold")) +
  ylab("Daily Returns/Daily VaR Forecasts") + 
  labs(caption = "Data source: Yahoo! Finance") +
  theme(plot.caption = element_text(face = "italic"))
  
# ************************************************************************************
# VaR Forecast on an Expanding Window, Forecasting Starting Date: 2018-12-31
# ************************************************************************************

garch_var       = ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(2,2)))
mdl_var         <- ugarchfit(spec = garch_var, data = SP500["/2018-12-31"]/100, solver = "hybrid")
std_res_var     <- residuals(mdl_var)/sigma(mdl_var)
g_var           <- gpd.fit(-std_res_var, threshold)
q_Z_var         <- threshold + g_var$mle[1] / g_var$mle[2] * (((1 - alpha) / g_var$rate) ^ (-g_var$mle[2]) - 1)
start_date      <- as.numeric(which(index(SP500) == "2018-12-31"))
VaR             <- as.vector(length(SP500["20181231/"]))

for (i in 1:(length(SP500["20181231/"])-1)) {
  data        <- SP500[(1):(start_date + i - 1)]/100
  garch       <- ugarchspec(variance.model = list(model = "eGARCH"))
  mdl         <- ugarchfit(garch, data, solver = "hybrid")
  forecast    <- ugarchforecast(mdl, n.ahead = 1)
  
  if (i %% 30 == 0) {
    data_gpd = SP500[1:(start_date + i)]
    garch_gpd     <- ugarchspec(variance.model = list(model = "eGARCH"))
    mdl_gpd       <- ugarchfit(garch_gpd, data_gpd, solver = "hybrid")
    gpd           <- gpd.fit(-residuals(mdl_gpd)/sigma(mdl_gpd), threshold)
    q_Z_var       <- threshold + gpd$mle[1] / gpd$mle[2] * (((1 - alpha) / gpd$rate) ^ (-gpd$mle[2]) - 1)
  }
  
  VaR[i] <- as.numeric(forecast@forecast$seriesFor) +  as.numeric(forecast@forecast$sigmaFor) * q_Z_var
}

VaR_forecast           <- data.frame(index(SP500["20190102/"]), VaR, data = SP500["20190102/"])
colnames(VaR_forecast) <- c("Date", "VaR", "SP500")

ggplot(VaR_forecast, aes(x = Date)) +
  geom_col(aes(y = SP500/100), width = 1, color = "lightslategray", fill = "lightslategray") +
  geom_line(aes(y = -VaR), color = "darkblue") +
  scale_y_continuous(labels = scales::percent) +
  xlab("") +
  theme_classic() +
  ggtitle("Daily VaR Forecast S&P 500", subtitle = "Expanding Window") +
  theme(plot.title = element_text(size = 18, face = "bold")) +
  ylab("Daily Returns/Daily VaR Forecasts") + 
  labs(caption = "Data source: Yahoo! Finance") +
  theme(plot.caption = element_text(face = "italic"))





