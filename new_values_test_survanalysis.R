library(tictoc)
library(parallel)
#A crude test to make sure that we don't alter the results of our simulation while optimising it
#Based on the results from Dani's original code

source('group dispersal version_commented_optimising.R')

set.seed(1) #would want this to change incrementally each run? - set to run number

ad.surv.intercept<-0.0002516054 
ad.surv.temp<-0.1702
ad.surv.pack<--0.1654
disp.intercept<-0.0064262
disp.pack.size<-0.10594
#IBI.average<-10.81
#IBI.temp<-1.02
pup.intercept<-0.87515
pup.pack.size<-0.04574
juv.surv.intercept<-18.61131
juv.surv.temp<--0.70572
juv.surv.litter<-0.5482
IBI.intercept<--17.25578
IBI.tempcoef<-0.91565
IBI.pups<-0.51976
#sets my mean temperature (change temp)
meantemp<-28.22864
#the actual mean temperature at present
actual.meantemp<-28.22864
number.of.packs<-9
n.steps = 100000
repeats = 16
burntime<-600
alphadieprob<-0.05149
variables<-16
#adsurvc<-c(0,repeats/16)
#adsurvtempc<-c(repeats/16,repeats/16*2)
#adsurvpackc<-c(repeats/16*2,repeats/16*3)
#dispinterceptc<-c(repeats/16*2,repeats/16*3)

#matrix of values
#a1<-c(ad.surv.intercept,ad.surv.temp,
#  ad.surv.pack,disp.intercept,disp.pack.size,IBI.intercept,IBI.tempcoef,IBI.pups,
#  pup.intercept,pup.pack.size,juv.surv.intercept,juv.surv.temp,juv.surv.litter,
#  meantemp,number.of.packs)

#m<-matrix(1,16,15)
#diag(m)<-1.01
#m[14,15]<-0.888888888888
#m[15,15]<-1.111111111
#m2<-t(replicate(16,a1))
#m3<-m2*m

cl <- makeCluster(no_cores)

tic()
results = run_simulation(ad.surv.intercept,ad.surv.temp,
                         ad.surv.pack,disp.intercept,disp.pack.size,IBI.intercept,IBI.tempcoef,IBI.pups,
                         pup.intercept,pup.pack.size,juv.surv.intercept,juv.surv.temp,juv.surv.litter,
                         meantemp,actual.meantemp,number.of.packs,n.steps,repeats, burntime)

toc()

View(results)

mean.pck = mean(results$mean.pck)
mean.pups = mean(results$mean.pups)
sd.pack = mean(results$sd.pack)
sd.IBI = mean(results$sd.IBI)
mean.IBI = mean(results$mean.IBI)
mean.blank = mean(results$mean.blank)
count.all.packs = mean(results$count.all.packs)



print("The following need to be TRUE when set.seed is set to 1")
print(all.equal(round(mean.pck,6),5.406149))
print(all.equal(round(mean.pups,6),3.645287))
print(all.equal(round(sd.pack,6),2.150649))
print(all.equal(round(sd.IBI,6),1.234829))
print(all.equal(round(mean.IBI,5),11.08547))
print(all.equal(round(mean.blank,6),7.222222))
print(all.equal(round(count.all.packs,6),8.945788))


