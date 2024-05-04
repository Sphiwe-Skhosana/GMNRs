###Importing and Pre-processing the data
data=read.csv("Data_Application1.csv")
xt=data$Urbanization ###non-parametric
zt=data$logEnergyUse_per_capita_Oil ##parametric
t=standardizer(xt,data$Code);data$logEnergyUse_per_capita_Oil=standardizer(zt,data$Code)
data=data[order(t),]
X=cbind(data$logEnergyUse_per_capita_Oil,sort(t))
x=X[,1];t=X[,2];
y=data$CO2_per_capita