###Importing and Pre-processing the data
data=read.csv("Data_Application2.csv")
xt=data$logGDP_per_capita ###non-parametric
zt=data$Renew_energy ##parametric
t=standardizer(xt,data$Code);data$Renew_energy=standardizer(zt,data$Code)
data=data[order(t),]
X=cbind(data$Renew_energy,sort(t))
x=X[,1];t=X[,2];
y=data$CO2_per_capita