# Time Series Analysis: ARIMA Modeling with Frequency Domain and Complex Root Decomposition
# Based on the report by Alexander Saines

# Load required libraries
library(forecast)
library(ggplot2)

# Helper functions for roots analysis
# Compute AR roots
arroots <- function(object) {
  if(!("Arima" %in% class(object)) & !("ar" %in% class(object)))
    stop("object must be of class Arima or ar")
  
  if("Arima" %in% class(object))
    parvec <- object$model$phi
  else
    parvec <- object$ar
  
  if(length(parvec) > 0) {
    last.nonzero <- max(which(abs(parvec) > 1e-08))
    if (last.nonzero > 0)
      return(structure(list(
        roots=polyroot(c(1,-parvec[1:last.nonzero])),
        type="AR"),
        class='armaroots'))
  }
  return(structure(list(roots=numeric(0), type="AR"),
                   class='armaroots'))
}

# Compute MA roots
maroots <- function(object) {
  if(!("Arima" %in% class(object)))
    stop("object must be of class Arima")
  
  parvec <- object$model$theta
  if(length(parvec) > 0) {
    last.nonzero <- max(which(abs(parvec) > 1e-08))
    if (last.nonzero > 0)
      return(structure(list(
        roots=polyroot(c(1,parvec[1:last.nonzero])),
        type="MA"),
        class='armaroots'))
  }
  return(structure(list(roots=numeric(0), type="MA"),
                   class='armaroots'))
}

# Plot ARMA roots
plot.armaroots <- function(x, xlab="Real", ylab="Imaginary",
                           main=paste("Inverse roots of", x$type,
                                     "characteristic polynomial"),
                           ...) {
  oldpar <- par(pty='s')
  on.exit(par(oldpar))
  plot(c(-1,1), c(-1,1), xlab=xlab, ylab=ylab,
       type="n", bty="n", xaxt="n", yaxt="n", main=main, ...)
  axis(1, at=c(-1,0,1), line=0.5, tck=-0.025)
  axis(2, at=c(-1,0,1), label=c("-i","0","i"),
       line=0.5, tck=-0.025)
  circx <- seq(-1,1,l=501)
  circy <- sqrt(1-circx^2)
  lines(c(circx,circx), c(circy,-circy), col='gray')
  lines(c(-2,2), c(0,0), col='gray')
  lines(c(0,0), c(-2,2), col='gray')
  if(length(x$roots) > 0) {
    inside <- abs(x$roots) > 1
    points(1/x$roots[inside], pch=19, col='black')
    if(sum(!inside) > 0)
      points(1/x$roots[!inside], pch=19, col='red')
  }
}

# Note: The arima.sim.freq2 function is assumed to be available
# If not available, you'll need to implement or import it

#--------------------------------------------------------------------------
# Exercise 1: Data Simulation using arima.sim.freq2
#--------------------------------------------------------------------------
# Given parameters
n <- 128
model <- list(
  Poles = c(0.9 * exp(2 * pi * 1i * 0.125), 0.8 * exp(2 * pi * 1i * 0.25)),
  Zeros = 0.7 * exp(2 * pi * 1i * 0.1875)
)
smodel <- list(
  Poles = 0.95,
  Zeros = 0.7 * exp(2 * pi * 1i * 0.1),
  S = 6
)
D <- 1
std <- 2
fs <- 4
unit <- "day"

# Simulate the data
result <- arima.sim.freq2(n, model, smodel, D, std = std, fs = fs, unit = unit)

# Display the plots
# Power Spectral Density (PSD) Plot
print(result$PSD)
print(result$ImpRes)
print(result$TS)
print(result$Configuration)

#--------------------------------------------------------------------------
# Exercise 2: Auto ARIMA Model
#--------------------------------------------------------------------------
# Convert to time series
ts1 <- ts(result$data$TimeSeries, frequency = 4)
autoplot(ts1)

## Fit auto.arima model
fitted_model <- auto.arima(ts1)
summary(fitted_model)

# Compare spectrum
ts1_d1 <- diff(ts1)

# Compute the spectrum using periodogram and autoregressive methods
s_period <- spectrum(ts1_d1, method = "pgram", fast = FALSE, log = "yes")
s_ar <- spectrum(ts1_d1, method = "ar", type = "b", log = "yes", 
                n.freq = (length(ts1_d1) / 2) + 1, order = 3)

# Create a data frame with the spectral densities
Spectral <- data.frame(
  freq = s_period$freq, 
  Period = s_period$spec, 
  AR = s_ar$spec[2:((length(ts1_d1) / 2) + 1)]
)

# Plot spectral densities
plot(ggplot(Spectral, aes(x=freq)) + 
     geom_line(aes(y = Period, colour = "periodogram")) + 
     geom_line(aes(y=AR, colour = "AR")))

# Compare angle
AR1 <- ar(ts1_d1, aic = FALSE, order = 3, method = c("ols"))
A1 <- c(1, -AR1$ar)
Fs = 4
AR_Model1 <- data.frame(
  Roots_fp = 1/polyroot(A1), 
  Radius = abs(1/polyroot(A1)), 
  Angle = Arg(1/polyroot(A1))/2/pi*Fs
)
print(AR_Model1)

#--------------------------------------------------------------------------
# Exercise 3: Manually Specified ARIMA Model
#--------------------------------------------------------------------------
# Fit the SARIMA model to the time series
SARIMA <- arima(ts1, order = c(4, 1, 2), 
               seasonal = list(order = c(1, 0, 2), period = 6), 
               method = "CSS")
summary(SARIMA)

# Analyze AR roots
A2 <- c(1, -SARIMA$coef[1:4])
AR_Model1 <- data.frame(
  Roots_fp = 1/polyroot(A2), 
  Radius = abs(1/polyroot(A2)), 
  Angle = Arg(1/polyroot(A2))/2/pi*Fs
)
print(AR_Model1)

# Analyze MA roots
B2 <- c(1, SARIMA$coef[5:6])
MA_Model1 <- data.frame(
  Roots_fp = 1/polyroot(B2), 
  Radius = abs(1/polyroot(B2)), 
  Angle = Arg(1/polyroot(B2))/2/pi*Fs
)
print(MA_Model1)

#--------------------------------------------------------------------------
# Exercise 4: Forecasting with the First 90 Samples
#--------------------------------------------------------------------------
# Split data into training and testing sets
train_data <- head(ts1, 90)
test_data <- tail(ts1, 38)

# First attempt - fit original structure
SARIMA2 <- arima(train_data, order = c(4, 1, 2), 
                seasonal = list(order = c(1, 0, 2), period = 6), 
                method = "CSS")
summary(SARIMA2)

# Plot inverse roots to check stationarity
par(mfrow=c(1,2))
plot(arroots(SARIMA2), main="Inverse AR roots")
plot(maroots(SARIMA2), main="Inverse MA roots")

# Structure without problems in roots
SARIMA2 <- arima(train_data, order = c(4, 1, 0), 
                seasonal = list(order = c(1, 0, 2), period = 6), 
                method = "CSS")
summary(SARIMA2)

# Plot inverse roots to check stationarity
par(mfrow=c(1,2))
plot(arroots(SARIMA2), main="Inverse AR roots")
plot(maroots(SARIMA2), main="Inverse MA roots")

# Forecast the future 38 samples using built-in forecast
forecasted_values <- forecast(SARIMA2, h = 38)
plot(autoplot(forecasted_values))

# Forecast using predictTS function (assuming it's available)
# If predictTS is not available, you'll need to implement or import it
# PredictTS(train_data, SARIMA2, 100, 38, 4, 0.05)

#--------------------------------------------------------------------------
# Exercise 5: Non-Seasonal ARIMA Model
#--------------------------------------------------------------------------
# Check spectrum of training data
ts1_d1 <- diff(train_data)

# Compute spectrum using periodogram and autoregressive methods (order 20 - Overfitting)
s_period <- spectrum(ts1_d1, method = "pgram", fast = FALSE, log = "yes")
s_ar <- spectrum(ts1_d1, method = "ar", type = "b", log = "yes", 
                n.freq = (length(ts1_d1) / 2) + 1, order = 20)

# Create data frame with spectral densities
Spectral <- data.frame(
  freq = s_period$freq, 
  Period = s_period$spec, 
  AR = s_ar$spec[2:((length(ts1_d1) / 2) + 1)]
)

# Plot spectral densities
plot(ggplot(Spectral, aes(x=freq)) + 
     geom_line(aes(y = Period, colour = "periodogram")) + 
     geom_line(aes(y=AR, colour = "AR")))

# First attempt - may be non-stationary
arima_train <- arima(train_data, order = c(4, 1, 2), method = "CSS")
par(mfrow=c(1,2))
plot(arroots(arima_train), main="Inverse AR roots")
plot(maroots(arima_train), main="Inverse MA roots")

# Check auto.arima on training dataset
fitted_model <- auto.arima(train_data)
summary(fitted_model)

# Fit ARIMA based on unit circle (stationary process), ground truth
arima_train2 <- arima(train_data, order = c(2, 1, 2), method = "CSS")
summary(arima_train2)

# Check stationarity
par(mfrow=c(1,2))
plot(arroots(arima_train2), main="Inverse AR roots")
plot(maroots(arima_train2), main="Inverse MA roots")

# Check the poles and zeros
A2 <- c(1, -arima_train2$coef[1:2])
AR_Model1 <- data.frame(
  Roots_fp = 1/polyroot(A2), 
  Radius = abs(1/polyroot(A2)), 
  Angle = Arg(1/polyroot(A2))/2/pi*Fs
)
print(AR_Model1)

B2 <- c(1, arima_train2$coef[3:4])
MA_Model1 <- data.frame(
  Roots_fp = 1/polyroot(B2), 
  Radius = abs(1/polyroot(B2)), 
  Angle = Arg(1/polyroot(B2))/2/pi*Fs
)
print(MA_Model1)

# Forecast using built-in forecast function
forecasted_values <- forecast(arima_train2, h = 38)
plot(autoplot(forecasted_values))

# Forecast using predictTS function (assuming it's available)
# PredictTS(train_data, arima_train2, 100, 38, 4, 0.05)
