##Imprting data and pre-processing
data=read.csv(data,paste0("CO2data-2014-Africa-new.csv"))
x=data$GDP.per.capita;x=standardizer(x)
y=data$CO2.emissions.per.capita;y=standardizer(y)
df=data.frame(x=x,y=y,country=data$Country.Name,code=data$Country.Code)[order(x),]
