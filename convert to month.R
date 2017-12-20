byday<-monthly.spreadsheet
byday.split <- split( byday , f = byday$ID )
View(byday.split[[1]])

aggregate(byday.split[[1]],by=month,FUN=mean)

a<-as.Date(byday$date, format = "%d/%m/%Y")
print(a)
byday$date<-a
byday$date<-format(byday$date,"%Y-%m-%d")
byday$date<-as.Date(byday$date, format = "%Y-%m-%d")
View(byday)

install.packages("xts")
library(xts)

install.packages("hydroTSM")
library(hydroTSM)
sumdat <- daily2monthly(byday.split, FUN=mean, na.rm=TRUE)
View(sumdat)

byday.split.zoo<-read.zoo(byday)

byday2 <- as.xts(byday)
sapply(byday,class)

death.m<-tapply(byday$death,paste(byday$ID,month(byday$date)),sum,na.rm=T)
View(death.m)
