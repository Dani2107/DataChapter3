library(tictoc)
library(parallel)
library(foreach)
library(doParallel)
library(plyr)
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
repeats = 1
burntime<-600
alphadieprob<-0.05149
variables<-16
#adsurvc<-c(0,repeats/16)
#adsurvtempc<-c(repeats/16,repeats/16*2)
#adsurvpackc<-c(repeats/16*2,repeats/16*3)
#dispinterceptc<-c(repeats/16*2,repeats/16*3)

#matrix of values
a1<-c(ad.surv.intercept,ad.surv.temp,
      ad.surv.pack,disp.intercept,disp.pack.size,IBI.intercept,IBI.tempcoef,IBI.pups,
      pup.intercept,pup.pack.size,juv.surv.intercept,juv.surv.temp,juv.surv.litter, 
      meantemp,actual.meantemp,number.of.packs,n.steps,repeats, burntime)
length(a1)

m<-matrix(1,16,19)
diag(m)<-1.01
m[16,16]<-1
m2<-t(replicate(16,a1))
m3<-as.matrix(m2*m)

View(m3)
m4<-as.data.frame(m3)


mean(m3[1,])

itx<-iter(m3, by='row')
View(itx)

tic()
results = run_simulation(ad.surv.intercept,ad.surv.temp,
                         ad.surv.pack,disp.intercept,disp.pack.size,IBI.intercept,IBI.tempcoef,IBI.pups,
                         pup.intercept,pup.pack.size,juv.surv.intercept,juv.surv.temp,juv.surv.litter,alphadieprob,
                         meantemp,actual.meantemp,number.of.packs,n.steps,repeats,burntime)

toc()

View(x)

View(results)

x <- m3
itx <- iter(x, by='row')
foreach(i=itx, .combine=c) %dopar% run_simulation(i[1])


cl <- makeCluster(12)
registerDoParallel(cl)

tic()
x <- foreach(i=1:16, .combine=rbind) %:%
     foreach(j=1:2, .packages=c("pracma","msm","matrixStats"),.combine=rbind) %dopar% run_simulation(ad.surv.intercept,ad.surv.temp,
                                                                                                         ad.surv.pack,disp.intercept,disp.pack.size,IBI.intercept,IBI.tempcoef,IBI.pups,
                                                                                                         pup.intercept,pup.pack.size,juv.surv.intercept,juv.surv.temp,juv.surv.litter,
                                                                                                         meantemp,actual.meantemp,number.of.packs,n.steps,repeats,burntime)

toc() 
warnings()


stopCluster(cl)

View(x)

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



df <-
  data.frame(
    ad.surv.intercept=rep(0.0002516054, 16),
    ad.surv.temp=0.1702,
    ad.surv.pack=-0.1654,
    disp.intercept=0.0064262,
    disp.pack.size=0.10594,
    IBI.intercept=-17.25578,
    IBI.tempcoef=0.91565,
    IBI.pups=0.51976,
    #IBI.average=10.81
    #IBI.temp=1.02
    pup.intercept=0.87515,
    pup.pack.size=0.04574,
    juv.surv.intercept=18.61131,
    juv.surv.temp=-0.70572,
    juv.surv.litter=0.5482,
    alphadieprob=0.05149,
    #sets my mean temperature (change temp)
    meantemp=28.22864,
    #the actual mean temperature at present
    actual.meantemp=28.22864,
    number.of.packs=9
  )

for(i in 1:(nrow(df)-1)){
  df[i, i] <- 1.01 * df[i, i]
}
#for(i in ((nrow(df)) + 1) : nrow(df)-1){
 # df[i, i - nrow(df)-1] <- 1.01 * df[i, i - nrow(df)-1]
#}



#for 2 iterations (guessing you want 1000 because it's stochastic)
its <- 1
df_big <- do.call(rbind, replicate(its, df, simplify = FALSE))


# Calculate the number of cores
no_cores <- detectCores()

# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, c("run_simulation","ad_surv","juv_surv","ad_disp","pup_no",
                    "IBI","alpha.die","all_death","adult_disp","lambda_func",
                    "actual.meantemp","number.of.packs","n.steps","repeats","burntime","df_big"))
clusterEvalQ(cl,library("pracma"))
clusterEvalQ(cl,library("msm"))
clusterEvalQ(cl,library("matrixStats"))                  

tic()
results <-
  parLapply(cl,1:nrow(df_big), function(x) 
        results = run_simulation(df_big[x, 1], df_big[x, 2],
                                 df_big[x, 3], df_big[x, 4],df_big[x, 5],df_big[x, 6],df_big[x, 7],df_big[x, 8],
                                 df_big[x, 9],df_big[x, 10],df_big[x, 11],df_big[x, 12],df_big[x, 13],df_big[x, 14],df_big[x, 15],
                                 actual.meantemp,number.of.packs,n.steps,repeats, burntime)
      
  )
  #lapply(1:nrow(df_big), function(x) 
   # results = run_simulation(df_big[x, 1], df_big[x, 2],
    #                         df_big[x, 3], df_big[x, 4],df_big[x, 5],df_big[x, 6],df_big[x, 7],df_big[x, 8],
     #                        df_big[x, 9],df_big[x, 10],df_big[x, 11],df_big[x, 12],df_big[x, 13],df_big[x, 14],df_big[x, 15],
      #                       actual.meantemp,number.of.packs,n.steps,repeats, burntime)
  #)  
toc()

r2<-do.call(rbind,results)

stopCluster(cl)

View(r2)
View(results)
results<-0
View(df)
results<-0
View(df_big)
df_big[,8]
