library(lme4)
library(nlme)
birth<-IBI_21Jul16
IBIglm<-glm(Month~maxtemp+pups3mo,data=birth,na.action=na.omit)
summary(IBIglm)
hist(birth$pups3mo)
NULLglm<-glm(time~1,data=birth)
anova(IBIglm,NULLglm)
