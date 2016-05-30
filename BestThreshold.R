
R <- c(24.12,35.81,10.85,11.31,26.56,27.90,30.15,33.92,34.47,36.47,37.31,41.12,43.10,45.46,47.26,52.26,53.50)

I <- c(4224,4438,4960.50,4271.75,4545.0,4747.00,5156.25,10722.0,4998.00,5818.50,7308.00,3750.00,5464.00,19990.0,10704.75,10308.0,15200.0)

Q <- data.frame(R,I)

library(ggplot2)

Q <- within(Q,Is<-smooth.spline(Q$R,Q$I,spar = 0.5)$y)

model <- lm(Q$I ~ Q$R + I(Q$R^2) + I(Q$R^3) + I(Q$R^4))

Q <- within(Q,Ip<-predict(model))

ggplot(Q)+geom_point(aes(R,I))+geom_line(aes(R,Is),col='blue')+geom_line(aes(R,Ip),col='red')
