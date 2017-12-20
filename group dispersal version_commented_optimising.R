#Dummy individual based model

library(pracma) 
library(profvis)
library(msm)
library(parallel)

# demographic functions
#function that calculates survival probablility of adults
ad_surv <- function(temp,size,a,b,c){
  u = a+b*temp+c*size
  u = exp(u)
  u = u/(1+u)
  return(u)}

#function that calculates survival probablility of juveniles
juv_surv<-function(temp,pups,a,b,c){
  u = a+b*temp+c*pups
  u = exp(u)
  u = u/(1+u)
  u = nthroot(u,9)
  return(u)}

#function that calculates dispersal probablility
ad_disp<-function(size,a,b){
  u = a+b*size
  u = exp(u)
  u = u/(1+u)
  return(u)}

#function that calculates number of pups
pup_no<-function(size,a,b){
  u = a+b*size
  u=round(u)
  return(u)}

#functions in loop
#alpha death - calculates if the alpha dies or not
alpha.die<-function(a,b,c,d){
  c<-(runif(a,0,1)<b)
  d<-length(c[c==TRUE])
}

#sum adult death - calculates the total number of adults (excluding the alpha) that die
all_death<-function(a,b,c,d){
  c<-(runif(a,0,1)<b)
  d<-length(c[c==FALSE])
}

#sum adult disp - calculates the number of adults that disperse
adult_disp<-function(a,b,c,d){
  c<-(runif(a,0,1)<b)
  if(length(c[c==TRUE])>=1){
    d<-sum(round(rtnorm((length(c[c==TRUE])),3,1,lower=0)))
  } else {
    d<-0
  }
  
}

#calculate lambda
lambda_func<-function(a,b){
  a<-a[b]/a[b-1]
  a[is.infinite(a)]=0
  print(a)
  #mean(a,na.rm=TRUE)
}

n.steps<-600


#profvis({
  ptm <- proc.time()
  


  

  repeats<-100 #number of times to repeat the loop
  dfall<-array(0,c(repeats+1,n.steps)) #number of dogs in the population t each timestep across runs
  #vectros of pack size, mean litter size, IBI, n timesteps empty, sd of pack size, IBI and pup number
  pck.vec.1<-pck.vec.2<-pck.vec.3<-pck.vec.4<-pck.vec.5<-pck.vec.6<-pck.vec.7<-pck.vec.8<-pck.vec.9<-litter.vec.1<-litter.vec.2<-litter.vec.3<-litter.vec.4<-litter.vec.5<-litter.vec.6<-litter.vec.7<-litter.vec.8<-litter.vec.9<- nbreed.vec.1<-nbreed.vec.2<-nbreed.vec.3<-nbreed.vec.4<-nbreed.vec.5<-nbreed.vec.6<-nbreed.vec.7<-nbreed.vec.8<-nbreed.vec.9<- blank.vec.1<-blank.vec.2<-blank.vec.3<-blank.vec.4<-blank.vec.5<-blank.vec.6<-blank.vec.7<-blank.vec.8<-blank.vec.9<-pck.vec.sd.1<-pck.vec.sd.2<-pck.vec.sd.3<-pck.vec.sd.4<-pck.vec.sd.5<-pck.vec.sd.6<-pck.vec.sd.7<-pck.vec.sd.8<-pck.vec.sd.9<-nbreed.vec.sd.1<-nbreed.vec.sd.2<-nbreed.vec.sd.3<-nbreed.vec.sd.4<-nbreed.vec.sd.5<-nbreed.vec.sd.6<-nbreed.vec.sd.7<-nbreed.vec.sd.8<-nbreed.vec.sd.9<-pup.sd.1<-pup.sd.2<-pup.sd.3<-pup.sd.4<-pup.sd.5<-pup.sd.6<-pup.sd.7<-pup.sd.8<-pup.sd.9<-rep(0,repeats)
  alpha.no.1<-alpha.no.2<-alpha.no.3<-alpha.no.4<-alpha.no.5<-alpha.no.6<-alpha.no.7<-alpha.no.8<-alpha.no.9<-1 #number of alphas in a particular matrix (pack)
  mean.pck<-sd.pack<-mean.pups<-sd.pups<-mean.IBI<-sd.IBI<-mean.blank<-sd.blank<-end.all.vec<-count.all.packs<-rep(0,repeats)
  dispdf<-array(0,c(18,n.steps+1))

  
  for(j in 1:repeats){
    # number of timesteps to run the model for ech time 
    n.steps <- 600
    # want to store the matrices, and the vector of population sizes
    ns.1 <- ns.2 <- ns.3 <- ns.4 <- ns.5 <- ns.6 <- ns.7 <- ns.8 <- ns.9 <- array(0,c(11,n.steps+1)) 
    # need to start with a population structure at time i=1
    ns.1[10,1]<-rpois(1,3)
    ns.2[10,1]<-rpois(1,3)
    ns.3[10,1]<-rpois(1,3)
    ns.4[10,1]<-rpois(1,3)
    ns.5[10,1]<-rpois(1,3)
    ns.6[10,1]<-rpois(1,3)
    ns.7[10,1]<-rpois(1,3)
    ns.8[10,1]<-rpois(1,3)
    ns.9[10,1]<-rpois(1,3)
    ns.1[11,1]<- ns.2[11,1]<- ns.3[11,1]<- ns.4[11,1]<- ns.5[11,1]<- ns.6[11,1]<- ns.7[11,1]<- ns.8[11,1]<- ns.9[11,1]<-1 
    #sets my litter size 1st run
    meanlitter<-round(rnorm(1,4,0.75))
    #sets my mean temperature (change temp)
    meantemp<-27.85
    #the actual mean temperature at present
    actual.meantemp<-27.85
    #sets when breeding will occurr first round
    nextbreed.1<-nextbreed.2<-nextbreed.3<-nextbreed.4<-nextbreed.5<-nextbreed.6<-nextbreed.7<-nextbreed.8<-nextbreed.9<-3
    k<-0 #counts from 0 each run
    k.vec <- temp.vec<-rep(0,n.steps) #stores k and temp in a vector
    #vector for pup number, pack number, IBI, number of packs, when there are no dogs left and the total number of dogs
    pups.vec.1<-pups.vec.2<-pups.vec.3<-pups.vec.4<-pups.vec.5<-pups.vec.6<-pups.vec.7<-pups.vec.8<-pups.vec.9<-pack.vec.1<-pack.vec.2<-pack.vec.3<-pack.vec.4<-pack.vec.5<-pack.vec.6<-pack.vec.7<-pack.vec.8<-pack.vec.9<-nextbreed.vec.1<-nextbreed.vec.2<-nextbreed.vec.3<-nextbreed.vec.4<-nextbreed.vec.5<-nextbreed.vec.6<-nextbreed.vec.7<-nextbreed.vec.8<-nextbreed.vec.9<-count.packs<-end.vec<-all.ad.packs<-rep(0,n.steps) 
    disp.ad.vec<-disp.ad.prev<-rep(0,9) #vector of dispersers and vector dispersers move into after 1 timestep
    #disp.ad.prev2<-rep(0,9) in case want to stay for 3 months
    #pack dispersal, the final run if all dog die, steps where there are no dogs and number of dogs in the pack
    packdisp.1<-packdisp.2<-packdisp.3<-packdisp.4<-packdisp.5<-packdisp.6<-packdisp.7<-packdisp.8<-packdisp.9<-last.1<-last.2<-last.3<-last.4<-last.5<-last.6<-last.7<-last.8<-last.9<-blank.1<-blank.2<-blank.3<-blank.4<-blank.5<-blank.6<-blank.7<-blank.8<-blank.9<-packnumb.1<-packnumb.2<-packnumb.3<-packnumb.4<-packnumb.5<-packnumb.6<-packnumb.7<-packnumb.8<-packnumb.9<-0 


      for (i in 1:n.steps){
      k<-k+1
      k.vec[i]<-k
      packnumb.1<-sum(ns.1[1:11,k-1])
      if(packnumb.1==0){
        blank.1<-blank.1+1
      }
      packnumb.2<-sum(ns.2[1:11,k-1])
      if(packnumb.2==0){
        blank.2<-blank.2+1
      }
      packnumb.3<-sum(ns.3[1:11,k-1])
      if(packnumb.3==0){
        blank.3<-blank.3+1
      }
      packnumb.4<-sum(ns.4[1:11,k-1])
      if(packnumb.4==0){
        blank.4<-blank.4+1
      }
      packnumb.5<-sum(ns.5[1:11,k-1])
      if(packnumb.5==0){
        blank.5<-blank.5+1
      }
      packnumb.6<-sum(ns.6[1:11,k-1])
      if(packnumb.6==0){
        blank.6<-blank.6+1
      }
      packnumb.7<-sum(ns.7[1:11,k-1])
      if(packnumb.7==0){
        blank.7<-blank.7+1
      }
      packnumb.8<-sum(ns.8[1:11,k-1])
      if(packnumb.8==0){
        blank.8<-blank.8+1
      }
      packnumb.9<-sum(ns.9[1:11,k-1])
      if(packnumb.9==0){
        blank.9<-blank.9+1
      }
      all.of.packs<-sum((c(packnumb.1,packnumb.2,packnumb.3,packnumb.4,packnumb.5,packnumb.6,packnumb.7,packnumb.8,packnumb.9)))
      all.ad.packs[i]<-all.of.packs
      if(all.of.packs==0){
        end.vec[i]<-k
      }
      
      count.packs[i]<-sum((c(packnumb.1,packnumb.2,packnumb.3,packnumb.4,packnumb.5,packnumb.6,packnumb.7,packnumb.8,packnumb.9))>0)
      #death
      death.1<-ns.1[11,k] #if there is an alpha death.n = 1, otherwise 0
      death.2<-ns.2[11,k]
      death.3<-ns.3[11,k]
      death.4<-ns.4[11,k]
      death.5<-ns.5[11,k]
      death.6<-ns.6[11,k]
      death.7<-ns.7[11,k]
      death.8<-ns.8[11,k]
      death.9<-ns.9[11,k]
      #define the size of the pack
      pack.size.1 <- sum(ns.1[1:11,i])
      pack.size.2 <- sum(ns.2[1:11,i])
      pack.size.3 <- sum(ns.3[1:11,i])
      pack.size.4 <- sum(ns.4[1:11,i])
      pack.size.5 <- sum(ns.5[1:11,i])
      pack.size.6 <- sum(ns.6[1:11,i])
      pack.size.7 <- sum(ns.7[1:11,i])
      pack.size.8 <- sum(ns.8[1:11,i])
      pack.size.9 <- sum(ns.9[1:11,i])
      pack.vec.1[i]<-pack.size.1
      pack.vec.2[i]<-pack.size.2
      pack.vec.3[i]<-pack.size.3
      pack.vec.4[i]<-pack.size.4
      pack.vec.5[i]<-pack.size.5
      pack.vec.6[i]<-pack.size.6
      pack.vec.7[i]<-pack.size.7
      pack.vec.8[i]<-pack.size.8
      pack.vec.9[i]<-pack.size.9
      #calculates litter size first run
      litter.size.1 <- pups.1 <- pup_no(pack.size.1,2.61734,0.19932) 
      litter.size.2 <- pups.2 <- pup_no(pack.size.2,2.61734,0.19932)
      litter.size.3 <- pups.3 <- pup_no(pack.size.3,2.61734,0.19932)
      litter.size.4 <- pups.4 <- pup_no(pack.size.4,2.61734,0.19932)
      litter.size.5 <- pups.5 <- pup_no(pack.size.5,2.61734,0.19932)
      litter.size.6 <- pups.6 <- pup_no(pack.size.6,2.61734,0.19932)
      litter.size.7 <- pups.7 <- pup_no(pack.size.7,2.61734,0.19932)
      litter.size.8 <- pups.8 <- pup_no(pack.size.8,2.61734,0.19932)
      litter.size.9 <- pups.9 <- pup_no(pack.size.9,2.61734,0.19932)
      pups.vec.1[i]<-pups.1
      pups.vec.2[i]<-pups.2
      pups.vec.3[i]<-pups.3
      pups.vec.4[i]<-pups.4
      pups.vec.5[i]<-pups.5
      pups.vec.6[i]<-pups.6
      pups.vec.7[i]<-pups.7
      pups.vec.8[i]<-pups.8
      pups.vec.9[i]<-pups.9
      # generate temperature
      temperature <- rnorm(1,meantemp,1)
      temp.vec[i]<-temperature
      #calculates the temperature over the previous 3 months - relevant for pup number d IBI
      # as gives temperature during the breeding period
      breedtemp<-((temp.vec[i]+temp.vec[i-1]+temp.vec[i-2])/3)
      #dispersers
      disp.all<-c(disp.ad.vec,disp.ad.prev) #disp.ad.prev.2 can add in if want dispersers to stay for 3 months
      dispdf[,i]<-disp.all
      #disp.ad.prev2<-disp.ad.prev
      disp.ad.prev<-disp.ad.vec
      dispsum<-sum(disp.all)
      dispfrac<-1/dispsum
      dispprob<-dispfrac*disp.all #probability of dispersal
      # adult survival
      surv.age.a.1 <- ad_surv(temperature,pack.size.1,9.47219,-0.22538,0.10292) 
      surv.age.a.2 <- ad_surv(temperature,pack.size.2,9.47219,-0.22538,0.10292)
      surv.age.a.3 <- ad_surv(temperature,pack.size.3,9.47219,-0.22538,0.10292)
      surv.age.a.4 <- ad_surv(temperature,pack.size.4,9.47219,-0.22538,0.10292)
      surv.age.a.5 <- ad_surv(temperature,pack.size.5,9.47219,-0.22538,0.10292)
      surv.age.a.6 <- ad_surv(temperature,pack.size.6,9.47219,-0.22538,0.10292)
      surv.age.a.7 <- ad_surv(temperature,pack.size.7,9.47219,-0.22538,0.10292)
      surv.age.a.8 <- ad_surv(temperature,pack.size.8,9.47219,-0.22538,0.10292)
      surv.age.a.9 <- ad_surv(temperature,pack.size.9,9.47219,-0.22538,0.10292)
      #number of non alpha adults
      ad.numb.1<-ns.1[10,k]
      ad.numb.2<-ns.2[10,k]
      ad.numb.3<-ns.3[10,k]
      ad.numb.4<-ns.4[10,k]
      ad.numb.5<-ns.5[10,k]
      ad.numb.6<-ns.6[10,k]
      ad.numb.7<-ns.7[10,k]
      ad.numb.8<-ns.8[10,k]
      ad.numb.9<-ns.9[10,k]
      #number of adults that die
      death.adult.1<-all_death(ad.numb.1,surv.age.a.1,d.ad.1,death.adult.1)
      death.adult.2<-all_death(ad.numb.2,surv.age.a.2,d.ad.2,death.adult.2)
      death.adult.3<-all_death(ad.numb.3,surv.age.a.3,d.ad.3,death.adult.3)
      death.adult.4<-all_death(ad.numb.4,surv.age.a.4,d.ad.4,death.adult.4)
      death.adult.5<-all_death(ad.numb.5,surv.age.a.5,d.ad.5,death.adult.5)
      death.adult.6<-all_death(ad.numb.6,surv.age.a.6,d.ad.6,death.adult.6)
      death.adult.7<-all_death(ad.numb.7,surv.age.a.7,d.ad.7,death.adult.7)
      death.adult.8<-all_death(ad.numb.8,surv.age.a.8,d.ad.8,death.adult.8)
      death.adult.9<-all_death(ad.numb.9,surv.age.a.9,d.ad.9,death.adult.9)
      #new number of adults after death
      ns.1[10,k+1]<-(ns.1[10,k]-death.adult.1)
      ns.2[10,k+1]<-(ns.2[10,k]-death.adult.2)
      ns.3[10,k+1]<-(ns.3[10,k]-death.adult.3)
      ns.4[10,k+1]<-(ns.4[10,k]-death.adult.4)
      ns.5[10,k+1]<-(ns.5[10,k]-death.adult.5)
      ns.6[10,k+1]<-(ns.6[10,k]-death.adult.6)
      ns.7[10,k+1]<-(ns.7[10,k]-death.adult.7)
      ns.8[10,k+1]<-(ns.8[10,k]-death.adult.8)
      ns.9[10,k+1]<-(ns.9[10,k]-death.adult.9)
      #dispersal
      disp.1<-ad_disp(pack.size.1,-4.49154,0.14) 
      disp.2<-ad_disp(pack.size.2,-4.49154,0.14)
      disp.3<-ad_disp(pack.size.3,-4.49154,0.14)
      disp.4<-ad_disp(pack.size.4,-4.49154,0.14)
      disp.5<-ad_disp(pack.size.5,-4.49154,0.14)
      disp.6<-ad_disp(pack.size.6,-4.49154,0.14)
      disp.7<-ad_disp(pack.size.7,-4.49154,0.14)
      disp.8<-ad_disp(pack.size.8,-4.49154,0.14)
      disp.9<-ad_disp(pack.size.9,-4.49154,0.14)
      #updates ad.numb to post death figures
      ad.numb.1<-ns.1[10,k+1]
      ad.numb.2<-ns.2[10,k+1]
      ad.numb.3<-ns.3[10,k+1] 
      ad.numb.4<-ns.4[10,k+1]
      ad.numb.5<-ns.5[10,k+1]
      ad.numb.6<-ns.6[10,k+1]
      ad.numb.7<-ns.7[10,k+1]
      ad.numb.8<-ns.8[10,k+1]
      ad.numb.9<-ns.9[10,k+1]
      #calculates the number of dispersing adults
      disp.adult.1<-adult_disp(ad.numb.1,disp.1,d.ad.1,disp.adult.1)
      if(disp.adult.1>ad.numb.1){
        disp.adult.1<-ad.numb.1
      }
      disp.adult.2<-adult_disp(ad.numb.2,disp.2,d.ad.2,disp.adult.2)
      if(disp.adult.2>ad.numb.2){
        disp.adult.2<-ad.numb.2
      }
      disp.adult.3<-adult_disp(ad.numb.3,disp.3,d.ad.3,disp.adult.3)
      if(disp.adult.3>ad.numb.3){
        disp.adult.3<-ad.numb.3
      }
      disp.adult.4<-adult_disp(ad.numb.4,disp.4,d.ad.4,disp.adult.4)
      if(disp.adult.4>ad.numb.4){
        disp.adult.4<-ad.numb.4
      }
      disp.adult.5<-adult_disp(ad.numb.5,disp.5,d.ad.5,disp.adult.5)
      if(disp.adult.5>ad.numb.5){
        disp.adult.5<-ad.numb.5
      }
      disp.adult.6<-adult_disp(ad.numb.6,disp.6,d.ad.6,disp.adult.6)
      if(disp.adult.6>ad.numb.6){
        disp.adult.6<-ad.numb.6
      }
      disp.adult.7<-adult_disp(ad.numb.7,disp.7,d.ad.7,disp.adult.7)
      if(disp.adult.7>ad.numb.7){
        disp.adult.7<-ad.numb.7
      }
      disp.adult.8<-adult_disp(ad.numb.8,disp.8,d.ad.8,disp.adult.8)
      if(disp.adult.8>ad.numb.8){
        disp.adult.8<-ad.numb.8
      }
      disp.adult.9<-adult_disp(ad.numb.9,disp.9,d.ad.9,disp.adult.9)
      if(disp.adult.9>ad.numb.9){
        disp.adult.9<-ad.numb.9
      }
      #adds dispersers to the dispersal vector
      disp.ad.vec[1]<-disp.adult.1+packdisp.1
      disp.ad.vec[2]<-disp.adult.2+packdisp.2
      disp.ad.vec[3]<-disp.adult.3+packdisp.3
      disp.ad.vec[4]<-disp.adult.4+packdisp.4
      disp.ad.vec[5]<-disp.adult.5+packdisp.5
      disp.ad.vec[6]<-disp.adult.6+packdisp.6
      disp.ad.vec[7]<-disp.adult.7+packdisp.7
      disp.ad.vec[8]<-disp.adult.8+packdisp.8
      disp.ad.vec[9]<-disp.adult.9+packdisp.9
      #object that gets filled with any pack dispersal numbers resulting from alpha death
      packdisp.1<-0
      packdisp.2<-0
      packdisp.3<-0
      packdisp.4<-0
      packdisp.5<-0
      packdisp.6<-0
      packdisp.7<-0
      packdisp.8<-0
      packdisp.9<-0
      #updates matrix with number of adults after dispersal
      ns.1[10,k+1]<-(ns.1[10,k+1]-disp.adult.1)
      ns.2[10,k+1]<-(ns.2[10,k+1]-disp.adult.2)
      ns.3[10,k+1]<-(ns.3[10,k+1]-disp.adult.3)
      ns.4[10,k+1]<-(ns.4[10,k+1]-disp.adult.4)
      ns.5[10,k+1]<-(ns.5[10,k+1]-disp.adult.5)
      ns.6[10,k+1]<-(ns.6[10,k+1]-disp.adult.6)
      ns.7[10,k+1]<-(ns.7[10,k+1]-disp.adult.7)
      ns.8[10,k+1]<-(ns.8[10,k+1]-disp.adult.8)
      ns.9[10,k+1]<-(ns.9[10,k+1]-disp.adult.9)
      #set juvenile survival
      surv.juv.1 <- juv_surv(temperature,litter.size.1,9.305655,-0.35286,0.5482)
      surv.juv.2 <- juv_surv(temperature,litter.size.2,9.305655,-0.35286,0.5482)
      surv.juv.3 <- juv_surv(temperature,litter.size.3,9.305655,-0.35286,0.5482)
      surv.juv.4 <- juv_surv(temperature,litter.size.4,9.305655,-0.35286,0.5482)
      surv.juv.5 <- juv_surv(temperature,litter.size.5,9.305655,-0.35286,0.5482)
      surv.juv.6 <- juv_surv(temperature,litter.size.6,9.305655,-0.35286,0.5482)
      surv.juv.7 <- juv_surv(temperature,litter.size.7,9.305655,-0.35286,0.5482)
      surv.juv.8 <- juv_surv(temperature,litter.size.8,9.305655,-0.35286,0.5482)
      surv.juv.9 <- juv_surv(temperature,litter.size.9,9.305655,-0.35286,0.5482)
      #gives which age class has juveniles in it
      nonzeroes.juv.1<-(which(ns.1[0:9,k]>=1))
      nonzeroes.juv.2<-(which(ns.2[0:9,k]>=1))
      nonzeroes.juv.3<-(which(ns.3[0:9,k]>=1))
      nonzeroes.juv.4<-(which(ns.4[0:9,k]>=1))
      nonzeroes.juv.5<-(which(ns.5[0:9,k]>=1))
      nonzeroes.juv.6<-(which(ns.6[0:9,k]>=1))
      nonzeroes.juv.7<-(which(ns.7[0:9,k]>=1))
      nonzeroes.juv.8<-(which(ns.8[0:9,k]>=1))
      nonzeroes.juv.9<-(which(ns.9[0:9,k]>=1))
      #number of juvenile age classes with dogs in
      length.juv.1<-length(nonzeroes.juv.1)
      length.juv.2<-length(nonzeroes.juv.2)
      length.juv.3<-length(nonzeroes.juv.3)
      length.juv.4<-length(nonzeroes.juv.4)
      length.juv.5<-length(nonzeroes.juv.5)
      length.juv.6<-length(nonzeroes.juv.6)
      length.juv.7<-length(nonzeroes.juv.7)
      length.juv.8<-length(nonzeroes.juv.8)
      length.juv.9<-length(nonzeroes.juv.9)
      #this calculates number of dogs that die each timestep. 
      if(length.juv.1 == 1) {
        nonzeroes.juv.first.1<-(nonzeroes.juv.1[1])
        if(nonzeroes.juv.first.1 < 9) {
          nonzeroes.juv.n.first.1<-(ns.1[nonzeroes.juv.first.1,k])
          death.juv.first.1<-all_death(nonzeroes.juv.n.first.1,surv.juv.1,c,d)
          ns.1[nonzeroes.juv.first.1+1,k+1]<-(ns.1[nonzeroes.juv.first.1,k]-death.juv.first.1) 
        } else {
          nonzeroes.juv.n.first.1<-(ns.1[nonzeroes.juv.first.1,k])
          death.juv.first.1<-all_death(nonzeroes.juv.n.first.1,surv.juv.1,c,d)
          ns.1[nonzeroes.juv.first.1+1,k+1]<-((ns.1[nonzeroes.juv.first.1,k]-death.juv.first.1) + ns.1[10,k])
        }
      }
      if(length.juv.1 > 1) {
        nonzeroes.juv.second.1<-(nonzeroes.juv.1[2])
        if (nonzeroes.juv.second.1 < 9) {
          nonzeroes.juv.n.second.1<-(ns.1[nonzeroes.juv.second.1,k])
          death.juv.second.1<-all_death(nonzeroes.juv.n.second.1,surv.juv.1,c,d)
          ns.1[nonzeroes.juv.second.1+1,k+1]<-(ns.1[nonzeroes.juv.second.1,k]-death.juv.second.1)
          nonzeroes.juv.first.1<-(nonzeroes.juv.1[1])
          nonzeroes.juv.n.first.1<-(ns.1[nonzeroes.juv.first.1,k])
          death.juv.first.1<-all_death(nonzeroes.juv.n.first.1,surv.juv.1,c,d)
          ns.1[nonzeroes.juv.first.1+1,k+1]<-(ns.1[nonzeroes.juv.first.1,k]-death.juv.first.1) 
        } else {
          nonzeroes.juv.n.second.1<-(ns.1[nonzeroes.juv.second.1,k])
          death.juv.second.1<-all_death(nonzeroes.juv.n.second.1,surv.juv.1,c,d)
          ns.1[nonzeroes.juv.second.1+1,k+1]<-((ns.1[nonzeroes.juv.second.1,k]-death.juv.second.1) + ns.1[10,k+1])
          nonzeroes.juv.first.1<-(nonzeroes.juv.1[1])
          nonzeroes.juv.n.first.1<-(ns.1[nonzeroes.juv.first.1,k])
          death.juv.first.1<-all_death(nonzeroes.juv.n.first.1,surv.juv.1,c,d)
          ns.1[nonzeroes.juv.first.1+1,k+1]<-(ns.1[nonzeroes.juv.first.1,k]-death.juv.first.1) 
        }
      }
      
      if(length.juv.2 == 1) {
        nonzeroes.juv.first.2<-(nonzeroes.juv.2[1])
        if(nonzeroes.juv.first.2 < 9) {
          nonzeroes.juv.n.first.2<-(ns.2[nonzeroes.juv.first.2,k])
          death.juv.first.2<-all_death(nonzeroes.juv.n.first.2,surv.juv.2,c,d)
          ns.2[nonzeroes.juv.first.2+1,k+1]<-(ns.2[nonzeroes.juv.first.2,k]-death.juv.first.2) 
        } else {
          nonzeroes.juv.n.first.2<-(ns.2[nonzeroes.juv.first.2,k])
          death.juv.first.2<-all_death(nonzeroes.juv.n.first.2,surv.juv.2,c,d)
          ns.2[nonzeroes.juv.first.2+1,k+1]<-((ns.2[nonzeroes.juv.first.2,k]-death.juv.first.2) + ns.2[10,k])
        }
      }
      if(length.juv.2 > 1) {
        nonzeroes.juv.second.2<-(nonzeroes.juv.2[2])
        if (nonzeroes.juv.second.2 < 9) {
          nonzeroes.juv.n.second.2<-(ns.2[nonzeroes.juv.second.2,k])
          death.juv.second.2<-all_death(nonzeroes.juv.n.second.2,surv.juv.2,c,d)
          ns.2[nonzeroes.juv.second.2+1,k+1]<-(ns.2[nonzeroes.juv.second.2,k]-death.juv.second.2)
          nonzeroes.juv.first.2<-(nonzeroes.juv.2[1])
          nonzeroes.juv.n.first.2<-(ns.2[nonzeroes.juv.first.2,k])
          death.juv.first.2<-all_death(nonzeroes.juv.n.first.2,surv.juv.2,c,d)
          ns.2[nonzeroes.juv.first.2+1,k+1]<-(ns.2[nonzeroes.juv.first.2,k]-death.juv.first.2) 
        } else {
          nonzeroes.juv.n.second.2<-(ns.2[nonzeroes.juv.second.2,k])
          death.juv.second.2<-all_death(nonzeroes.juv.n.second.2,surv.juv.2,c,d)
          ns.2[nonzeroes.juv.second.2+1,k+1]<-((ns.2[nonzeroes.juv.second.2,k]-death.juv.second.2) + ns.2[10,k+1])
          nonzeroes.juv.first.2<-(nonzeroes.juv.2[1])
          nonzeroes.juv.n.first.2<-(ns.2[nonzeroes.juv.first.2,k])
          death.juv.first.2<-all_death(nonzeroes.juv.n.first.2,surv.juv.2,c,d)
          ns.2[nonzeroes.juv.first.2+1,k+1]<-(ns.2[nonzeroes.juv.first.2,k]-death.juv.first.2) 
        }
      }
      
      if(length.juv.3 == 1) {
        nonzeroes.juv.first.3<-(nonzeroes.juv.3[1])
        if(nonzeroes.juv.first.3 < 9) {
          nonzeroes.juv.n.first.3<-(ns.3[nonzeroes.juv.first.3,k])
          death.juv.first.3<-all_death(nonzeroes.juv.n.first.3,surv.juv.3,c,d)
          ns.3[nonzeroes.juv.first.3+1,k+1]<-(ns.3[nonzeroes.juv.first.3,k]-death.juv.first.3) 
        } else {
          nonzeroes.juv.n.first.3<-(ns.3[nonzeroes.juv.first.3,k])
          death.juv.first.3<-all_death(nonzeroes.juv.n.first.3,surv.juv.3,c,d)
          ns.3[nonzeroes.juv.first.3+1,k+1]<-((ns.3[nonzeroes.juv.first.3,k]-death.juv.first.3) + ns.3[10,k])
        }
      }
      if(length.juv.3 > 1) {
        nonzeroes.juv.second.3<-(nonzeroes.juv.3[2])
        if (nonzeroes.juv.second.3 < 9) {
          nonzeroes.juv.n.second.3<-(ns.3[nonzeroes.juv.second.3,k])
          death.juv.second.3<-all_death(nonzeroes.juv.n.second.3,surv.juv.3,c,d)
          ns.3[nonzeroes.juv.second.3+1,k+1]<-(ns.3[nonzeroes.juv.second.3,k]-death.juv.second.3)
          nonzeroes.juv.first.3<-(nonzeroes.juv.3[1])
          nonzeroes.juv.n.first.3<-(ns.3[nonzeroes.juv.first.3,k])
          death.juv.first.3<-all_death(nonzeroes.juv.n.first.3,surv.juv.3,c,d)
          ns.3[nonzeroes.juv.first.3+1,k+1]<-(ns.3[nonzeroes.juv.first.3,k]-death.juv.first.3) 
        } else {
          nonzeroes.juv.n.second.3<-(ns.3[nonzeroes.juv.second.3,k])
          death.juv.second.3<-all_death(nonzeroes.juv.n.second.3,surv.juv.3,c,d)
          ns.3[nonzeroes.juv.second.3+1,k+1]<-((ns.3[nonzeroes.juv.second.3,k]-death.juv.second.3) + ns.3[10,k+1])
          nonzeroes.juv.first.3<-(nonzeroes.juv.3[1])
          nonzeroes.juv.n.first.3<-(ns.3[nonzeroes.juv.first.3,k])
          death.juv.first.3<-all_death(nonzeroes.juv.n.first.3,surv.juv.3,c,d)
          ns.3[nonzeroes.juv.first.3+1,k+1]<-(ns.3[nonzeroes.juv.first.3,k]-death.juv.first.3) 
        }
      }
      
      if(length.juv.4 == 1) {
        nonzeroes.juv.first.4<-(nonzeroes.juv.4[1])
        if(nonzeroes.juv.first.4 < 9) {
          nonzeroes.juv.n.first.4<-(ns.4[nonzeroes.juv.first.4,k])
          death.juv.first.4<-all_death(nonzeroes.juv.n.first.4,surv.juv.4,c,d)
          ns.4[nonzeroes.juv.first.4+1,k+1]<-(ns.4[nonzeroes.juv.first.4,k]-death.juv.first.4) 
        } else {
          nonzeroes.juv.n.first.4<-(ns.4[nonzeroes.juv.first.4,k])
          death.juv.first.4<-all_death(nonzeroes.juv.n.first.4,surv.juv.4,c,d)
          ns.4[nonzeroes.juv.first.4+1,k+1]<-((ns.4[nonzeroes.juv.first.4,k]-death.juv.first.4) + ns.4[10,k])
        }
      }
      if(length.juv.4 > 1) {
        nonzeroes.juv.second.4<-(nonzeroes.juv.4[2])
        if (nonzeroes.juv.second.4 < 9) {
          nonzeroes.juv.n.second.4<-(ns.4[nonzeroes.juv.second.4,k])
          death.juv.second.4<-all_death(nonzeroes.juv.n.second.4,surv.juv.4,c,d)
          ns.4[nonzeroes.juv.second.4+1,k+1]<-(ns.4[nonzeroes.juv.second.4,k]-death.juv.second.4)
          nonzeroes.juv.first.4<-(nonzeroes.juv.4[1])
          nonzeroes.juv.n.first.4<-(ns.4[nonzeroes.juv.first.4,k])
          death.juv.first.4<-all_death(nonzeroes.juv.n.first.4,surv.juv.4,c,d)
          ns.4[nonzeroes.juv.first.4+1,k+1]<-(ns.4[nonzeroes.juv.first.4,k]-death.juv.first.4) 
        } else {
          nonzeroes.juv.n.second.4<-(ns.4[nonzeroes.juv.second.4,k])
          death.juv.second.4<-all_death(nonzeroes.juv.n.second.4,surv.juv.4,c,d)
          ns.4[nonzeroes.juv.second.4+1,k+1]<-((ns.4[nonzeroes.juv.second.4,k]-death.juv.second.4) + ns.4[10,k+1])
          nonzeroes.juv.first.4<-(nonzeroes.juv.4[1])
          nonzeroes.juv.n.first.4<-(ns.4[nonzeroes.juv.first.4,k])
          death.juv.first.4<-all_death(nonzeroes.juv.n.first.4,surv.juv.4,c,d)
          ns.4[nonzeroes.juv.first.4+1,k+1]<-(ns.4[nonzeroes.juv.first.4,k]-death.juv.first.4) 
        }
      }
      
      if(length.juv.5 == 1) {
        nonzeroes.juv.first.5<-(nonzeroes.juv.5[1])
        if(nonzeroes.juv.first.5 < 9) {
          nonzeroes.juv.n.first.5<-(ns.5[nonzeroes.juv.first.5,k])
          death.juv.first.5<-all_death(nonzeroes.juv.n.first.5,surv.juv.5,c,d)
          ns.5[nonzeroes.juv.first.5+1,k+1]<-(ns.5[nonzeroes.juv.first.5,k]-death.juv.first.5) 
        } else {
          nonzeroes.juv.n.first.5<-(ns.5[nonzeroes.juv.first.5,k])
          death.juv.first.5<-all_death(nonzeroes.juv.n.first.5,surv.juv.5,c,d)
          ns.5[nonzeroes.juv.first.5+1,k+1]<-((ns.5[nonzeroes.juv.first.5,k]-death.juv.first.5) + ns.5[10,k])
        }
      }
      if(length.juv.5 > 1) {
        nonzeroes.juv.second.5<-(nonzeroes.juv.5[2])
        if (nonzeroes.juv.second.5 < 9) {
          nonzeroes.juv.n.second.5<-(ns.5[nonzeroes.juv.second.5,k])
          death.juv.second.5<-all_death(nonzeroes.juv.n.second.5,surv.juv.5,c,d)
          ns.5[nonzeroes.juv.second.5+1,k+1]<-(ns.5[nonzeroes.juv.second.5,k]-death.juv.second.5)
          nonzeroes.juv.first.5<-(nonzeroes.juv.5[1])
          nonzeroes.juv.n.first.5<-(ns.5[nonzeroes.juv.first.5,k])
          death.juv.first.5<-all_death(nonzeroes.juv.n.first.5,surv.juv.5,c,d)
          ns.5[nonzeroes.juv.first.5+1,k+1]<-(ns.5[nonzeroes.juv.first.5,k]-death.juv.first.5) 
        } else {
          nonzeroes.juv.n.second.5<-(ns.5[nonzeroes.juv.second.5,k])
          death.juv.second.5<-all_death(nonzeroes.juv.n.second.5,surv.juv.5,c,d)
          ns.5[nonzeroes.juv.second.5+1,k+1]<-((ns.5[nonzeroes.juv.second.5,k]-death.juv.second.5) + ns.5[10,k+1])
          nonzeroes.juv.first.5<-(nonzeroes.juv.5[1])
          nonzeroes.juv.n.first.5<-(ns.5[nonzeroes.juv.first.5,k])
          death.juv.first.5<-all_death(nonzeroes.juv.n.first.5,surv.juv.5,c,d)
          ns.5[nonzeroes.juv.first.5+1,k+1]<-(ns.5[nonzeroes.juv.first.5,k]-death.juv.first.5) 
        }
      }
      
      if(length.juv.6 == 1) {
        nonzeroes.juv.first.6<-(nonzeroes.juv.6[1])
        if(nonzeroes.juv.first.6 < 9) {
          nonzeroes.juv.n.first.6<-(ns.6[nonzeroes.juv.first.6,k])
          death.juv.first.6<-all_death(nonzeroes.juv.n.first.6,surv.juv.6,c,d)
          ns.6[nonzeroes.juv.first.6+1,k+1]<-(ns.6[nonzeroes.juv.first.6,k]-death.juv.first.6) 
        } else {
          nonzeroes.juv.n.first.6<-(ns.6[nonzeroes.juv.first.6,k])
          death.juv.first.6<-all_death(nonzeroes.juv.n.first.6,surv.juv.6,c,d)
          ns.6[nonzeroes.juv.first.6+1,k+1]<-((ns.6[nonzeroes.juv.first.6,k]-death.juv.first.6) + ns.6[10,k])
        }
      }
      if(length.juv.6 > 1) {
        nonzeroes.juv.second.6<-(nonzeroes.juv.6[2])
        if (nonzeroes.juv.second.6 < 9) {
          nonzeroes.juv.n.second.6<-(ns.6[nonzeroes.juv.second.6,k])
          death.juv.second.6<-all_death(nonzeroes.juv.n.second.6,surv.juv.6,c,d)
          ns.6[nonzeroes.juv.second.6+1,k+1]<-(ns.6[nonzeroes.juv.second.6,k]-death.juv.second.6)
          nonzeroes.juv.first.6<-(nonzeroes.juv.6[1])
          nonzeroes.juv.n.first.6<-(ns.6[nonzeroes.juv.first.6,k])
          death.juv.first.6<-all_death(nonzeroes.juv.n.first.6,surv.juv.6,c,d)
          ns.6[nonzeroes.juv.first.6+1,k+1]<-(ns.6[nonzeroes.juv.first.6,k]-death.juv.first.6) 
        } else {
          nonzeroes.juv.n.second.6<-(ns.6[nonzeroes.juv.second.6,k])
          death.juv.second.6<-all_death(nonzeroes.juv.n.second.6,surv.juv.6,c,d)
          ns.6[nonzeroes.juv.second.6+1,k+1]<-((ns.6[nonzeroes.juv.second.6,k]-death.juv.second.6) + ns.6[10,k+1])
          nonzeroes.juv.first.6<-(nonzeroes.juv.6[1])
          nonzeroes.juv.n.first.6<-(ns.6[nonzeroes.juv.first.6,k])
          death.juv.first.6<-all_death(nonzeroes.juv.n.first.6,surv.juv.6,c,d)
          ns.6[nonzeroes.juv.first.6+1,k+1]<-(ns.6[nonzeroes.juv.first.6,k]-death.juv.first.6) 
        }
      }
      
      if(length.juv.7 == 1) {
        nonzeroes.juv.first.7<-(nonzeroes.juv.7[1])
        if(nonzeroes.juv.first.7 < 9) {
          nonzeroes.juv.n.first.7<-(ns.7[nonzeroes.juv.first.7,k])
          death.juv.first.7<-all_death(nonzeroes.juv.n.first.7,surv.juv.7,c,d)
          ns.7[nonzeroes.juv.first.7+1,k+1]<-(ns.7[nonzeroes.juv.first.7,k]-death.juv.first.7) 
        } else {
          nonzeroes.juv.n.first.7<-(ns.7[nonzeroes.juv.first.7,k])
          death.juv.first.7<-all_death(nonzeroes.juv.n.first.7,surv.juv.7,c,d)
          ns.7[nonzeroes.juv.first.7+1,k+1]<-((ns.7[nonzeroes.juv.first.7,k]-death.juv.first.7) + ns.7[10,k])
        }
      }
      if(length.juv.7 > 1) {
        nonzeroes.juv.second.7<-(nonzeroes.juv.7[2])
        if (nonzeroes.juv.second.7 < 9) {
          nonzeroes.juv.n.second.7<-(ns.7[nonzeroes.juv.second.7,k])
          death.juv.second.7<-all_death(nonzeroes.juv.n.second.7,surv.juv.7,c,d)
          ns.7[nonzeroes.juv.second.7+1,k+1]<-(ns.7[nonzeroes.juv.second.7,k]-death.juv.second.7)
          nonzeroes.juv.first.7<-(nonzeroes.juv.7[1])
          nonzeroes.juv.n.first.7<-(ns.7[nonzeroes.juv.first.7,k])
          death.juv.first.7<-all_death(nonzeroes.juv.n.first.7,surv.juv.7,c,d)
          ns.7[nonzeroes.juv.first.7+1,k+1]<-(ns.7[nonzeroes.juv.first.7,k]-death.juv.first.7) 
        } else {
          nonzeroes.juv.n.second.7<-(ns.7[nonzeroes.juv.second.7,k])
          death.juv.second.7<-all_death(nonzeroes.juv.n.second.7,surv.juv.7,c,d)
          ns.7[nonzeroes.juv.second.7+1,k+1]<-((ns.7[nonzeroes.juv.second.7,k]-death.juv.second.7) + ns.7[10,k+1])
          nonzeroes.juv.first.7<-(nonzeroes.juv.7[1])
          nonzeroes.juv.n.first.7<-(ns.7[nonzeroes.juv.first.7,k])
          death.juv.first.7<-all_death(nonzeroes.juv.n.first.7,surv.juv.7,c,d)
          ns.7[nonzeroes.juv.first.7+1,k+1]<-(ns.7[nonzeroes.juv.first.7,k]-death.juv.first.7) 
        }
      }
      
      if(length.juv.8 == 1) {
        nonzeroes.juv.first.8<-(nonzeroes.juv.8[1])
        if(nonzeroes.juv.first.8 < 9) {
          nonzeroes.juv.n.first.8<-(ns.8[nonzeroes.juv.first.8,k])
          death.juv.first.8<-all_death(nonzeroes.juv.n.first.8,surv.juv.8,c,d)
          ns.8[nonzeroes.juv.first.8+1,k+1]<-(ns.8[nonzeroes.juv.first.8,k]-death.juv.first.8) 
        } else {
          nonzeroes.juv.n.first.8<-(ns.8[nonzeroes.juv.first.8,k])
          death.juv.first.8<-all_death(nonzeroes.juv.n.first.8,surv.juv.8,c,d)
          ns.8[nonzeroes.juv.first.8+1,k+1]<-((ns.8[nonzeroes.juv.first.8,k]-death.juv.first.8) + ns.8[10,k])
        }
      }
      if(length.juv.8 > 1) {
        nonzeroes.juv.second.8<-(nonzeroes.juv.8[2])
        if (nonzeroes.juv.second.8 < 9) {
          nonzeroes.juv.n.second.8<-(ns.8[nonzeroes.juv.second.8,k])
          death.juv.second.8<-all_death(nonzeroes.juv.n.second.8,surv.juv.8,c,d)
          ns.8[nonzeroes.juv.second.8+1,k+1]<-(ns.8[nonzeroes.juv.second.8,k]-death.juv.second.8)
          nonzeroes.juv.first.8<-(nonzeroes.juv.8[1])
          nonzeroes.juv.n.first.8<-(ns.8[nonzeroes.juv.first.8,k])
          death.juv.first.8<-all_death(nonzeroes.juv.n.first.8,surv.juv.8,c,d)
          ns.8[nonzeroes.juv.first.8+1,k+1]<-(ns.8[nonzeroes.juv.first.8,k]-death.juv.first.8) 
        } else {
          nonzeroes.juv.n.second.8<-(ns.8[nonzeroes.juv.second.8,k])
          death.juv.second.8<-all_death(nonzeroes.juv.n.second.8,surv.juv.8,c,d)
          ns.8[nonzeroes.juv.second.8+1,k+1]<-((ns.8[nonzeroes.juv.second.8,k]-death.juv.second.8) + ns.8[10,k+1])
          nonzeroes.juv.first.8<-(nonzeroes.juv.8[1])
          nonzeroes.juv.n.first.8<-(ns.8[nonzeroes.juv.first.8,k])
          death.juv.first.8<-all_death(nonzeroes.juv.n.first.8,surv.juv.8,c,d)
          ns.8[nonzeroes.juv.first.8+1,k+1]<-(ns.8[nonzeroes.juv.first.8,k]-death.juv.first.8) 
        }
      }
      
      if(length.juv.9 == 1) {
        nonzeroes.juv.first.9<-(nonzeroes.juv.9[1])
        if(nonzeroes.juv.first.9 < 9) {
          nonzeroes.juv.n.first.9<-(ns.9[nonzeroes.juv.first.9,k])
          death.juv.first.9<-all_death(nonzeroes.juv.n.first.9,surv.juv.9,c,d)
          ns.9[nonzeroes.juv.first.9+1,k+1]<-(ns.9[nonzeroes.juv.first.9,k]-death.juv.first.9) 
        } else {
          nonzeroes.juv.n.first.9<-(ns.9[nonzeroes.juv.first.9,k])
          death.juv.first.9<-all_death(nonzeroes.juv.n.first.9,surv.juv.9,c,d)
          ns.9[nonzeroes.juv.first.9+1,k+1]<-((ns.9[nonzeroes.juv.first.9,k]-death.juv.first.9) + ns.9[10,k])
        }
      }
      if(length.juv.9 > 1) {
        nonzeroes.juv.second.9<-(nonzeroes.juv.9[2])
        if (nonzeroes.juv.second.9 < 9) {
          nonzeroes.juv.n.second.9<-(ns.9[nonzeroes.juv.second.9,k])
          death.juv.second.9<-all_death(nonzeroes.juv.n.second.9,surv.juv.9,c,d)
          ns.9[nonzeroes.juv.second.9+1,k+1]<-(ns.9[nonzeroes.juv.second.9,k]-death.juv.second.9)
          nonzeroes.juv.first.9<-(nonzeroes.juv.9[1])
          nonzeroes.juv.n.first.9<-(ns.9[nonzeroes.juv.first.9,k])
          death.juv.first.9<-all_death(nonzeroes.juv.n.first.9,surv.juv.9,c,d)
          ns.9[nonzeroes.juv.first.9+1,k+1]<-(ns.9[nonzeroes.juv.first.9,k]-death.juv.first.9) 
        } else {
          nonzeroes.juv.n.second.9<-(ns.9[nonzeroes.juv.second.9,k])
          death.juv.second.9<-all_death(nonzeroes.juv.n.second.9,surv.juv.9,c,d)
          ns.9[nonzeroes.juv.second.9+1,k+1]<-((ns.9[nonzeroes.juv.second.9,k]-death.juv.second.9) + ns.9[10,k+1])
          nonzeroes.juv.first.9<-(nonzeroes.juv.9[1])
          nonzeroes.juv.n.first.9<-(ns.9[nonzeroes.juv.first.9,k])
          death.juv.first.9<-all_death(nonzeroes.juv.n.first.9,surv.juv.9,c,d)
          ns.9[nonzeroes.juv.first.9+1,k+1]<-(ns.9[nonzeroes.juv.first.9,k]-death.juv.first.9) 
        }
      }
      # now do the probability of recruitment - ie do they give birth or not, varies with temp
      if (i == nextbreed.1) {
        pb.1 <- 1
      } else {
        pb.1 <- 0
      }
      if (i == nextbreed.2) {
        pb.2 <- 1
      } else {
        pb.2 <- 0
      }
      if (i == nextbreed.3) {
        pb.3 <- 1
      } else {
        pb.3 <- 0
      }
      if (i == nextbreed.4) {
        pb.4 <- 1
      } else {
        pb.4 <- 0
      }
      if (i == nextbreed.5) {
        pb.5 <- 1
      } else {
        pb.5 <- 0
      }
      if (i == nextbreed.6) {
        pb.6 <- 1
      } else {
        pb.6 <- 0
      }
      if (i == nextbreed.7) {
        pb.7 <- 1
      } else {
        pb.7 <- 0
      }
      if (i == nextbreed.8) {
        pb.8 <- 1
      } else {
        pb.8 <- 0
      }
      if (i == nextbreed.9) {
        pb.9 <- 1
      } else {
        pb.9 <- 0
      }
      # Pick a the next step (i) when breeding occurs.
      if (pb.1 == 1) {
        tempdiff <- breedtemp - actual.meantemp
        litterdiff.1 <- pups.1 - meanlitter
        nextbreed.sum.1<-nextbreed.1
        nextbreed.1 <- round(rnorm(1, mean = (10.89 + (tempdiff*0.8209) + (litterdiff.1/2)), 1) + i) 
        nextbreed.vec.1[i]<-nextbreed.1-nextbreed.sum.1
      }
      if (pb.2 == 1) {
        tempdiff <- breedtemp - actual.meantemp
        litterdiff.2 <- pups.2 - meanlitter
        nextbreed.sum.2<-nextbreed.2
        nextbreed.2 <- round(rnorm(1, mean = (10.89 + (tempdiff*0.8209) + (litterdiff.2/2)), 1) + i)
        nextbreed.vec.2[i]<-nextbreed.2-nextbreed.sum.2
      }
      if (pb.3 == 1) {
        tempdiff <- breedtemp - actual.meantemp
        litterdiff.3 <- pups.3 - meanlitter
        nextbreed.sum.3<-nextbreed.3
        nextbreed.3 <- round(rnorm(1, mean = (10.89 + (tempdiff*0.8209) + (litterdiff.3/2)), 1) + i)
        nextbreed.vec.3[i]<-nextbreed.3-nextbreed.sum.3
      }
      if (pb.4 == 1) {
        tempdiff <- breedtemp - actual.meantemp
        litterdiff.4 <- pups.4 - meanlitter
        nextbreed.sum.4<-nextbreed.4
        nextbreed.4 <- round(rnorm(1, mean = (10.89 + (tempdiff*0.8209) + (litterdiff.4/2)), 1) + i)
        nextbreed.vec.4[i]<-nextbreed.4-nextbreed.sum.4
      }
      if (pb.5 == 1) {
        tempdiff <- breedtemp - actual.meantemp
        litterdiff.5 <- pups.5 - meanlitter
        nextbreed.sum.5<-nextbreed.5
        nextbreed.5 <- round(rnorm(1, mean = (10.89 + (tempdiff*0.8209) + (litterdiff.5/2)), 1) + i)
        nextbreed.vec.5[i]<-nextbreed.5-nextbreed.sum.5
      }
      if (pb.6 == 1) {
        tempdiff <- breedtemp - actual.meantemp
        litterdiff.6 <- pups.6 - meanlitter
        nextbreed.sum.6<-nextbreed.6
        nextbreed.6 <- round(rnorm(1, mean = (10.89 + (tempdiff*0.8209) + (litterdiff.6/2)), 1) + i)
        nextbreed.vec.6[i]<-nextbreed.6-nextbreed.sum.6
      }
      if (pb.7 == 1) {
        tempdiff <- breedtemp - actual.meantemp
        litterdiff.7 <- pups.7 - meanlitter
        nextbreed.sum.7<-nextbreed.7
        nextbreed.7 <- round(rnorm(1, mean = (10.89 + (tempdiff*0.8209) + (litterdiff.7/2)), 1) + i)
        nextbreed.vec.7[i]<-nextbreed.7-nextbreed.sum.7
      }
      if (pb.8 == 1) {
        tempdiff <- breedtemp - actual.meantemp
        litterdiff.8 <- pups.8 - meanlitter
        nextbreed.sum.8<-nextbreed.8
        nextbreed.8 <- round(rnorm(1, mean = (10.89 + (tempdiff*0.8209) + (litterdiff.8/2)), 1) + i)
        nextbreed.vec.8[i]<-nextbreed.8-nextbreed.sum.8
      }
      if (pb.9 == 1) {
        tempdiff <- breedtemp - actual.meantemp
        litterdiff.9 <- pups.9 - meanlitter
        nextbreed.sum.9<-nextbreed.9
        nextbreed.9 <- round(rnorm(1, mean = (10.89 + (tempdiff*0.8209) + (litterdiff.9/2)), 1) + i)
        nextbreed.vec.9[i]<-nextbreed.9-nextbreed.sum.9
      }
      # now number of pups produced
      pups.1 <- pup_no(pack.size.1,2.7453,0.1797) 
      litter.size.1<-pb.1*pups.1
      if(ns.1[11,k]==1){
        ns.1[1,k+1] <- pb.1*pups.1
      } else {
        ns.1[1,k+1] <- 0
      }
      pups.2 <- pup_no(pack.size.2,2.7453,0.1797)
      litter.size.2<-pb.2*pups.2
      if(ns.2[11,k]==1){
        ns.2[1,k+1] <- pb.2*pups.2
      } else {
        ns.2[1,k+1] <- 0
      }
      pups.3 <- pup_no(pack.size.3,2.7453,0.1797)
      litter.size.3<-pb.3*pups.3
      if(ns.3[11,k]==1){
        ns.3[1,k+1] <- pb.3*pups.3
      } else {
        ns.3[1,k+1] <- 0
      }
      pups.4 <- pup_no(pack.size.4,2.7453,0.1797) 
      litter.size.4<-pb.4*pups.4
      if(ns.4[11,k]==1){
        ns.4[1,k+1] <- pb.4*pups.4
      } else {
        ns.4[1,k+1] <- 0
      }
      pups.5 <- pup_no(pack.size.5,2.7453,0.1797)
      litter.size.5<-pb.5*pups.5
      if(ns.5[11,k]==1){
        ns.5[1,k+1] <- pb.5*pups.5
      } else {
        ns.5[1,k+1] <- 0
      }
      pups.6 <- pup_no(pack.size.6,2.7453,0.1797)
      litter.size.6<-pb.6*pups.6
      if(ns.6[11,k]==1){
        ns.6[1,k+1] <- pb.6*pups.6
      } else {
        ns.6[1,k+1] <- 0
      }
      pups.7 <- pup_no(pack.size.7,2.7453,0.1797) 
      litter.size.7<-pb.7*pups.7
      if(ns.7[11,k]==1){
        ns.7[1,k+1] <- pb.7*pups.7
      } else {
        ns.7[1,k+1] <- 0
      }
      pups.8 <- pup_no(pack.size.8,2.61734,0.1797)
      litter.size.8<-pb.8*pups.8
      if(ns.8[11,k]==1){
        ns.8[1,k+1] <- pb.8*pups.8
      } else {
        ns.8[1,k+1] <- 0
      }
      pups.9 <- pup_no(pack.size.9,2.7453,0.1797)
      litter.size.9<-pb.9*pups.9
      if(ns.9[11,k]==1){
        ns.9[1,k+1] <- pb.9*pups.9
      } else {
        ns.9[1,k+1] <- 0
      }
      #alpha death dynamic
      # now dominant suvival
      #this sets whether the alpha has died (0) or not (1) and assigns it to 'death'
      if (alpha.no.1==1){
        alpha.no.1<-death.1<-alpha.die(1,surv.age.a.1,c,d)
      }
      ns.1[11,k+1] <- alpha.no.1
      
      if (alpha.no.2==1){
        alpha.no.2<-death.2<-alpha.die(1,surv.age.a.2,c,d)
      }
      ns.2[11,k+1] <- alpha.no.2
      
      if (alpha.no.3==1){
        alpha.no.3<-death.3<-alpha.die(1,surv.age.a.3,c,d)
      }
      ns.3[11,k+1] <- alpha.no.3
      
      if (alpha.no.4==1){
        alpha.no.4<-death.4<-alpha.die(1,surv.age.a.4,c,d)
      }
      ns.4[11,k+1] <- alpha.no.4
      
      if (alpha.no.5==1){
        alpha.no.5<-death.5<-alpha.die(1,surv.age.a.5,c,d)
      }
      ns.5[11,k+1] <- alpha.no.5
      
      if (alpha.no.6==1){
        alpha.no.6<-death.6<-alpha.die(1,surv.age.a.6,c,d)
      }
      ns.6[11,k+1] <- alpha.no.6
      
      if (alpha.no.7==1){
        alpha.no.7<-death.7<-alpha.die(1,surv.age.a.7,c,d)
      }
      ns.7[11,k+1] <- alpha.no.7
      
      if (alpha.no.8==1){
        alpha.no.8<-death.8<-alpha.die(1,surv.age.a.8,c,d)
      }
      ns.8[11,k+1] <- alpha.no.8
      
      if (alpha.no.9==1){
        alpha.no.9<-death.9<-alpha.die(1,surv.age.a.9,c,d)
      }
      ns.9[11,k+1] <- alpha.no.9
#this determines what will happen if the alpha has died. Either the whole pack disperses
# or another adult takes over. If there are juveniles present during whole pack dispersal they die
      if (death.1 == 0) {
        dth <- runif(1,0,100)
        if (dth > 40) {
          ns.1[11,k+1]<-0
          last.1<-k
          #if the whole pack has dispersed + the matrix is empty then this picks which dispersers 
          #fill the matrix (form a new pack). One becomes an alpha. if theres >1 the others are 
          #non alpha adults
          if (sum(disp.all) > 0) {
            replace.1<-sample(disp.all, replace=FALSE, prob=dispprob,1)
            pos.all.1<-which(replace.1==disp.all)
            pos.1<-sample(pos.all.1,1)
            disp.all[pos.1]<-0
            if (pos.1 <= 9) {
              disp.ad.vec[pos.1]<-0
            }
          } else {
            if (sum(disp.all) == 0) {
              replace.1<-0
            }  
          }
          if (replace.1 == 1) {
            ns.1[1:10,k+1]<-0
            ns.1[11,k+1]<-1
            death.1<-1
            alpha.no.1<-1
          } else {
            if (replace.1 > 1) {
              ns.1[1:9,k+1]<-0
              ns.1[10,k+1]<-replace.1-1
              ns.1[11,k+1]<-1
              death.1<-1
              alpha.no.1<-1
            } else {
              if (replace.1 == 0) {
                ns.1[1:11,k+1]<-0
                death.1<-1
                alpha.no.1<-0
              }
            }
          }
          packdisp.1<-disp.adult.1+sum(ns.1[9:10,i])
          
          
        } else {
          
          while (death.1 == 0) {
            if (ns.1[10, k] >= 1) {
              
              ns.1[10, k+1] <- ns.1[10, k] - 1
              death.1 <- 1
              alpha.no.1 <-1
              ns.1[11,k+1]<-death.1
              
            } else {
              
              # These two statements stop the loop going on forever.
              if (ns.1[10, k] == 0) {
                if (sum(disp.all) > 0) {
                  replace.1<-sample(disp.all, replace=FALSE, prob=dispprob,1)
                  pos.all.1<-which(replace.1==disp.all)
                  pos.1<-sample(pos.all.1,1)
                  disp.all[pos.1]<-0
                  if (pos.1 <= 9) {
                    disp.ad.vec[pos.1]<-0
                  }
                } else{
                  if (sum(disp.all) == 0) {
                    replace.1<-0
                  }   
                }
                if (replace.1 == 1) {
                  ns.1[1:10,k+1]<-0
                  ns.1[11,k+1]<-1
                  alpha.no.1 <- 1
                } else {
                  if (replace.1 > 1) {
                    ns.1[1:9,k+1]<-0
                    ns.1[10,k+1]<-replace.1-1
                    ns.1[11,k+1]<-1
                    death.1<-1
                    alpha.no.1<-1
                  } else {
                    if (replace.1 == 0) {
                      ns.1[1:11,k+1]<-0
                      death.1<-1
                      alpha.no.1<-0
                    }
                  }
                }
                last.1<-k
                
                
                
              }
            }
          }
        }
      }
      if (death.2 == 0) {
        dth <- runif(1,0,100)
        if (dth > 40) {
          ns.2[11,k+1]<-0
          last.2<-k
          if (sum(disp.all) > 0) {
            replace.2<-sample(disp.all, replace=FALSE, prob=dispprob,1)
            pos.all.2<-which(replace.2==disp.all)
            pos.2<-sample(pos.all.2,1)
            disp.all[pos.2]<-0
            if (pos.2 <= 9) {
              disp.ad.vec[pos.2]<-0
            }
          } else {
            if (sum(disp.all) == 0) {
              replace.2<-0
            }  
          }
          if (replace.2 == 1) {
            ns.2[1:10,k+1]<-0
            ns.2[11,k+1]<-1
            death.2<-1
            alpha.no.2<-1
          } else {
            if (replace.2 > 1) {
              ns.2[1:9,k+1]<-0
              ns.2[10,k+1]<-replace.2-1
              ns.2[11,k+1]<-1
              death.2<-1
              alpha.no.2<-1
            } else {
              if (replace.2 == 0) {
                ns.2[1:11,k+1]<-0
                death.2<-1
                alpha.no.2<-0
              }
            }
          }
          packdisp.2<-disp.adult.2+sum(ns.2[9:10,i])
          
          
        } else {
          
          while (death.2 == 0) {
            if (ns.2[10, k] >= 1) {
              
              ns.2[10, k+1] <- ns.2[10, k] - 1
              death.2 <- 1
              alpha.no.2 <-1
              ns.2[11,k+1]<-death.2
              
            } else {
              
              # These two statements stop the loop going on forever.
              if (ns.2[10, k] == 0) {
                if (sum(disp.all) > 0) {
                  replace.2<-sample(disp.all, replace=FALSE, prob=dispprob,1)
                  pos.all.2<-which(replace.2==disp.all)
                  pos.2<-sample(pos.all.2,1)
                  disp.all[pos.2]<-0
                  if (pos.2 <= 9) {
                    disp.ad.vec[pos.2]<-0
                  }
                } else{
                  if (sum(disp.all) == 0) {
                    replace.2<-0
                  }   
                }
                if (replace.2 == 1) {
                  ns.2[1:10,k+1]<-0
                  ns.2[11,k+1]<-1
                  alpha.no.2 <- 1
                } else {
                  if (replace.2 > 1) {
                    ns.2[1:9,k+1]<-0
                    ns.2[10,k+1]<-replace.2-1
                    ns.2[11,k+1]<-1
                    death.2<-1
                    alpha.no.2<-1
                  } else {
                    if (replace.2 == 0) {
                      ns.2[1:11,k+1]<-0
                      death.2<-1
                      alpha.no.2<-0
                    }
                  }
                }
                last.2<-k
                
                
                
              }
            }
          }
        }
      }
      if (death.3 == 0) {
        dth <- runif(1,0,100)
        if (dth > 40) {
          ns.3[11,k+1]<-0
          last.3<-k
          if (sum(disp.all) > 0) {
            replace.3<-sample(disp.all, replace=FALSE, prob=dispprob,1)
            pos.all.3<-which(replace.3==disp.all)
            pos.3<-sample(pos.all.3,1)
            disp.all[pos.3]<-0
            if (pos.3 <= 9) {
              disp.ad.vec[pos.3]<-0
            }
          } else {
            if (sum(disp.all) == 0) {
              replace.3<-0
            }  
          }
          if (replace.3 == 1) {
            ns.3[1:10,k+1]<-0
            ns.3[11,k+1]<-1
            death.3<-1
            alpha.no.3<-1
          } else {
            if (replace.3 > 1) {
              ns.3[1:9,k+1]<-0
              ns.3[10,k+1]<-replace.3-1
              ns.3[11,k+1]<-1
              death.3<-1
              alpha.no.3<-1
            } else {
              if (replace.3 == 0) {
                ns.3[1:11,k+1]<-0
                death.3<-1
                alpha.no.3<-0
              }
            }
          }
          packdisp.3<-disp.adult.3+sum(ns.3[9:10,i])
          
          
        } else {
          
          while (death.3 == 0) {
            if (ns.3[10, k] >= 1) {
              
              ns.3[10, k+1] <- ns.3[10, k] - 1
              death.3 <- 1
              alpha.no.3 <-1
              ns.3[11,k+1]<-death.3
              
            } else {
              
              # These two statements stop the loop going on forever.
              if (ns.3[10, k] == 0) {
                if (sum(disp.all) > 0) {
                  replace.3<-sample(disp.all, replace=FALSE, prob=dispprob,1)
                  pos.all.3<-which(replace.3==disp.all)
                  pos.3<-sample(pos.all.3,1)
                  disp.all[pos.3]<-0
                  if (pos.3 <= 9) {
                    disp.ad.vec[pos.3]<-0
                  }
                } else{
                  if (sum(disp.all) == 0) {
                    replace.3<-0
                  }   
                }
                if (replace.3 == 1) {
                  ns.3[1:10,k+1]<-0
                  ns.3[11,k+1]<-1
                  alpha.no.3 <- 1
                } else {
                  if (replace.3 > 1) {
                    ns.3[1:9,k+1]<-0
                    ns.3[10,k+1]<-replace.3-1
                    ns.3[11,k+1]<-1
                    death.3<-1
                    alpha.no.3<-1
                  } else {
                    if (replace.3 == 0) {
                      ns.3[1:11,k+1]<-0
                      death.3<-1
                      alpha.no.3<-0
                    }
                  }
                }
                last.3<-k
                
                
                
              }
            }
          }
        }
      }
      if (death.4 == 0) {
        dth <- runif(1,0,100)
        if (dth > 40) {
          ns.4[11,k+1]<-0
          last.4<-k
          if (sum(disp.all) > 0) {
            replace.4<-sample(disp.all, replace=FALSE, prob=dispprob,1)
            pos.all.4<-which(replace.4==disp.all)
            pos.4<-sample(pos.all.4,1)
            disp.all[pos.4]<-0
            if (pos.4 <= 9) {
              disp.ad.vec[pos.4]<-0
            }
          } else {
            if (sum(disp.all) == 0) {
              replace.4<-0
            }  
          }
          if (replace.4 == 1) {
            ns.4[1:10,k+1]<-0
            ns.4[11,k+1]<-1
            death.4<-1
            alpha.no.4<-1
          } else {
            if (replace.4 > 1) {
              ns.4[1:9,k+1]<-0
              ns.4[10,k+1]<-replace.4-1
              ns.4[11,k+1]<-1
              death.4<-1
              alpha.no.4<-1
            } else {
              if (replace.4 == 0) {
                ns.4[1:11,k+1]<-0
                death.4<-1
                alpha.no.4<-0
              }
            }
          }
          packdisp.4<-disp.adult.4+sum(ns.4[9:10,i])
          
          
        } else {
          
          while (death.4 == 0) {
            if (ns.4[10, k] >= 1) {
              
              ns.4[10, k+1] <- ns.4[10, k] - 1
              death.4 <- 1
              alpha.no.4 <-1
              ns.4[11,k+1]<-death.4
              
            } else {
              
              # These two statements stop the loop going on forever.
              if (ns.4[10, k] == 0) {
                if (sum(disp.all) > 0) {
                  replace.4<-sample(disp.all, replace=FALSE, prob=dispprob,1)
                  pos.all.4<-which(replace.4==disp.all)
                  pos.4<-sample(pos.all.4,1)
                  disp.all[pos.4]<-0
                  if (pos.4 <= 9) {
                    disp.ad.vec[pos.4]<-0
                  }
                } else{
                  if (sum(disp.all) == 0) {
                    replace.4<-0
                  }   
                }
                if (replace.4 == 1) {
                  ns.4[1:10,k+1]<-0
                  ns.4[11,k+1]<-1
                  alpha.no.4 <- 1
                } else {
                  if (replace.4 > 1) {
                    ns.4[1:9,k+1]<-0
                    ns.4[10,k+1]<-replace.4-1
                    ns.4[11,k+1]<-1
                    death.4<-1
                    alpha.no.4<-1
                  } else {
                    if (replace.4 == 0) {
                      ns.4[1:11,k+1]<-0
                      death.4<-1
                      alpha.no.4<-0
                    }
                  }
                }
                last.4<-k
                
                
                
              }
            }
          }
        }
      }
      if (death.5 == 0) {
        dth <- runif(1,0,100)
        if (dth > 40) {
          ns.5[11,k+1]<-0
          last.5<-k
          if (sum(disp.all) > 0) {
            replace.5<-sample(disp.all, replace=FALSE, prob=dispprob,1)
            pos.all.5<-which(replace.5==disp.all)
            pos.5<-sample(pos.all.5,1)
            disp.all[pos.5]<-0
            if (pos.5 <= 9) {
              disp.ad.vec[pos.5]<-0
            }
          } else {
            if (sum(disp.all) == 0) {
              replace.5<-0
            }  
          }
          if (replace.5 == 1) {
            ns.5[1:10,k+1]<-0
            ns.5[11,k+1]<-1
            death.5<-1
            alpha.no.5<-1
          } else {
            if (replace.5 > 1) {
              ns.5[1:9,k+1]<-0
              ns.5[10,k+1]<-replace.5-1
              ns.5[11,k+1]<-1
              death.5<-1
              alpha.no.5<-1
            } else {
              if (replace.5 == 0) {
                ns.5[1:11,k+1]<-0
                death.5<-1
                alpha.no.5<-0
              }
            }
          }
          packdisp.5<-disp.adult.5+sum(ns.5[9:10,i])
          
          
        } else {
          
          while (death.5 == 0) {
            if (ns.5[10, k] >= 1) {
              
              ns.5[10, k+1] <- ns.5[10, k] - 1
              death.5 <- 1
              alpha.no.5 <-1
              ns.5[11,k+1]<-death.5
              
            } else {
              
              # These two statements stop the loop going on forever.
              if (ns.5[10, k] == 0) {
                if (sum(disp.all) > 0) {
                  replace.5<-sample(disp.all, replace=FALSE, prob=dispprob,1)
                  pos.all.5<-which(replace.5==disp.all)
                  pos.5<-sample(pos.all.5,1)
                  disp.all[pos.5]<-0
                  if (pos.5 <= 9) {
                    disp.ad.vec[pos.5]<-0
                  }
                } else{
                  if (sum(disp.all) == 0) {
                    replace.5<-0
                  }   
                }
                if (replace.5 == 1) {
                  ns.5[1:10,k+1]<-0
                  ns.5[11,k+1]<-1
                  alpha.no.5 <- 1
                } else {
                  if (replace.5 > 1) {
                    ns.5[1:9,k+1]<-0
                    ns.5[10,k+1]<-replace.5-1
                    ns.5[11,k+1]<-1
                    death.5<-1
                    alpha.no.5<-1
                  } else {
                    if (replace.5 == 0) {
                      ns.5[1:11,k+1]<-0
                      death.5<-1
                      alpha.no.5<-0
                    }
                  }
                }
                last.5<-k
                
                
                
              }
            }
          }
        }
      }
      if (death.6 == 0) {
        dth <- runif(1,0,100)
        if (dth > 40) {
          ns.6[11,k+1]<-0
          last.6<-k
          if (sum(disp.all) > 0) {
            replace.6<-sample(disp.all, replace=FALSE, prob=dispprob,1)
            pos.all.6<-which(replace.6==disp.all)
            pos.6<-sample(pos.all.6,1)
            disp.all[pos.6]<-0
            if (pos.6 <= 9) {
              disp.ad.vec[pos.6]<-0
            }
          } else {
            if (sum(disp.all) == 0) {
              replace.6<-0
            }  
          }
          if (replace.6 == 1) {
            ns.6[1:10,k+1]<-0
            ns.6[11,k+1]<-1
            death.6<-1
            alpha.no.6<-1
          } else {
            if (replace.6 > 1) {
              ns.6[1:9,k+1]<-0
              ns.6[10,k+1]<-replace.6-1
              ns.6[11,k+1]<-1
              death.6<-1
              alpha.no.6<-1
            } else {
              if (replace.6 == 0) {
                ns.6[1:11,k+1]<-0
                death.6<-1
                alpha.no.6<-0
              }
            }
          }
          packdisp.6<-disp.adult.6+sum(ns.6[9:10,i])
          
          
        } else {
          
          while (death.6 == 0) {
            if (ns.6[10, k] >= 1) {
              
              ns.6[10, k+1] <- ns.6[10, k] - 1
              death.6 <- 1
              alpha.no.6 <-1
              ns.6[11,k+1]<-death.6
              
            } else {
              
              # These two statements stop the loop going on forever.
              if (ns.6[10, k] == 0) {
                if (sum(disp.all) > 0) {
                  replace.6<-sample(disp.all, replace=FALSE, prob=dispprob,1)
                  pos.all.6<-which(replace.6==disp.all)
                  pos.6<-sample(pos.all.6,1)
                  disp.all[pos.6]<-0
                  if (pos.6 <= 9) {
                    disp.ad.vec[pos.6]<-0
                  }
                } else{
                  if (sum(disp.all) == 0) {
                    replace.6<-0
                  }   
                }
                if (replace.6 == 1) {
                  ns.6[1:10,k+1]<-0
                  ns.6[11,k+1]<-1
                  alpha.no.6 <- 1
                } else {
                  if (replace.6 > 1) {
                    ns.6[1:9,k+1]<-0
                    ns.6[10,k+1]<-replace.6-1
                    ns.6[11,k+1]<-1
                    death.6<-1
                    alpha.no.6<-1
                  } else {
                    if (replace.6 == 0) {
                      ns.6[1:11,k+1]<-0
                      death.6<-1
                      alpha.no.6<-0
                    }
                  }
                }
                last.6<-k
                
                
                
              }
            }
          }
        }
      }
      if (death.7 == 0) {
        dth <- runif(1,0,100)
        if (dth > 40) {
          ns.7[11,k+1]<-0
          last.7<-k
          if (sum(disp.all) > 0) {
            replace.7<-sample(disp.all, replace=FALSE, prob=dispprob,1)
            pos.all.7<-which(replace.7==disp.all)
            pos.7<-sample(pos.all.7,1)
            disp.all[pos.7]<-0
            if (pos.7 <= 9) {
              disp.ad.vec[pos.7]<-0
            }
          } else {
            if (sum(disp.all) == 0) {
              replace.7<-0
            }  
          }
          if (replace.7 == 1) {
            ns.7[1:10,k+1]<-0
            ns.7[11,k+1]<-1
            death.7<-1
            alpha.no.7<-1
          } else {
            if (replace.7 > 1) {
              ns.7[1:9,k+1]<-0
              ns.7[10,k+1]<-replace.7-1
              ns.7[11,k+1]<-1
              death.7<-1
              alpha.no.7<-1
            } else {
              if (replace.7 == 0) {
                ns.7[1:11,k+1]<-0
                death.7<-1
                alpha.no.7<-0
              }
            }
          }
          packdisp.7<-disp.adult.7+sum(ns.7[9:10,i])
          
          
        } else {
          
          while (death.7 == 0) {
            if (ns.7[10, k] >= 1) {
              
              ns.7[10, k+1] <- ns.7[10, k] - 1
              death.7 <- 1
              alpha.no.7 <-1
              ns.7[11,k+1]<-death.7
              
            } else {
              
              # These two statements stop the loop going on forever.
              if (ns.7[10, k] == 0) {
                if (sum(disp.all) > 0) {
                  replace.7<-sample(disp.all, replace=FALSE, prob=dispprob,1)
                  pos.all.7<-which(replace.7==disp.all)
                  pos.7<-sample(pos.all.7,1)
                  disp.all[pos.7]<-0
                  if (pos.7 <= 9) {
                    disp.ad.vec[pos.7]<-0
                  }
                } else{
                  if (sum(disp.all) == 0) {
                    replace.7<-0
                  }   
                }
                if (replace.7 == 1) {
                  ns.7[1:10,k+1]<-0
                  ns.7[11,k+1]<-1
                  alpha.no.7 <- 1
                } else {
                  if (replace.7 > 1) {
                    ns.7[1:9,k+1]<-0
                    ns.7[10,k+1]<-replace.7-1
                    ns.7[11,k+1]<-1
                    death.7<-1
                    alpha.no.7<-1
                  } else {
                    if (replace.7 == 0) {
                      ns.7[1:11,k+1]<-0
                      death.7<-1
                      alpha.no.7<-0
                    }
                  }
                }
                last.7<-k
                
                
                
              }
            }
          }
        }
      }
      if (death.8 == 0) {
        dth <- runif(1,0,100)
        if (dth > 40) {
          ns.8[11,k+1]<-0
          last.8<-k
          if (sum(disp.all) > 0) {
            replace.8<-sample(disp.all, replace=FALSE, prob=dispprob,1)
            pos.all.8<-which(replace.8==disp.all)
            pos.8<-sample(pos.all.8,1)
            disp.all[pos.8]<-0
            if (pos.8 <= 9) {
              disp.ad.vec[pos.8]<-0
            }
          } else {
            if (sum(disp.all) == 0) {
              replace.8<-0
            }  
          }
          if (replace.8 == 1) {
            ns.8[1:10,k+1]<-0
            ns.8[11,k+1]<-1
            death.8<-1
            alpha.no.8<-1
          } else {
            if (replace.8 > 1) {
              ns.8[1:9,k+1]<-0
              ns.8[10,k+1]<-replace.8-1
              ns.8[11,k+1]<-1
              death.8<-1
              alpha.no.8<-1
            } else {
              if (replace.8 == 0) {
                ns.8[1:11,k+1]<-0
                death.8<-1
                alpha.no.8<-0
              }
            }
          }
          packdisp.8<-disp.adult.8+sum(ns.8[9:10,i])
          
          
        } else {
          
          while (death.8 == 0) {
            if (ns.8[10, k] >= 1) {
              
              ns.8[10, k+1] <- ns.8[10, k] - 1
              death.8 <- 1
              alpha.no.8 <-1
              ns.8[11,k+1]<-death.8
              
            } else {
              
              # These two statements stop the loop going on forever.
              if (ns.8[10, k] == 0) {
                if (sum(disp.all) > 0) {
                  replace.8<-sample(disp.all, replace=FALSE, prob=dispprob,1)
                  pos.all.8<-which(replace.8==disp.all)
                  pos.8<-sample(pos.all.8,1)
                  disp.all[pos.8]<-0
                  if (pos.8 <= 9) {
                    disp.ad.vec[pos.8]<-0
                  }
                } else{
                  if (sum(disp.all) == 0) {
                    replace.8<-0
                  }   
                }
                if (replace.8 == 1) {
                  ns.8[1:10,k+1]<-0
                  ns.8[11,k+1]<-1
                  alpha.no.8 <- 1
                } else {
                  if (replace.8 > 1) {
                    ns.8[1:9,k+1]<-0
                    ns.8[10,k+1]<-replace.8-1
                    ns.8[11,k+1]<-1
                    death.8<-1
                    alpha.no.8<-1
                  } else {
                    if (replace.8 == 0) {
                      ns.8[1:11,k+1]<-0
                      death.8<-1
                      alpha.no.8<-0
                    }
                  }
                }
                last.8<-k
                
                
                
              }
            }
          }
        }
      }
      if (death.9 == 0) {
        dth <- runif(1,0,100)
        if (dth > 40) {
          ns.9[11,k+1]<-0
          last.9<-k
          if (sum(disp.all) > 0) {
            replace.9<-sample(disp.all, replace=FALSE, prob=dispprob,1)
            pos.all.9<-which(replace.9==disp.all)
            pos.9<-sample(pos.all.9,1)
            disp.all[pos.9]<-0
            if (pos.9 <= 9) {
              disp.ad.vec[pos.9]<-0
            }
          } else {
            if (sum(disp.all) == 0) {
              replace.9<-0
            }  
          }
          if (replace.9 == 1) {
            ns.9[1:10,k+1]<-0
            ns.9[11,k+1]<-1
            death.9<-1
            alpha.no.9<-1
          } else {
            if (replace.9 > 1) {
              ns.9[1:9,k+1]<-0
              ns.9[10,k+1]<-replace.9-1
              ns.9[11,k+1]<-1
              death.9<-1
              alpha.no.9<-1
            } else {
              if (replace.9 == 0) {
                ns.9[1:11,k+1]<-0
                death.9<-1
                alpha.no.9<-0
              }
            }
          }
          packdisp.9<-disp.adult.9+sum(ns.9[9:10,i])
          
          
        } else {
          
          while (death.9 == 0) {
            if (ns.9[10, k] >= 1) {
              
              ns.9[10, k+1] <- ns.9[10, k] - 1
              death.9 <- 1
              alpha.no.9 <-1
              ns.9[11,k+1]<-death.9
              
            } else {
              
              # These two statements stop the loop going on forever.
              if (ns.9[10, k] == 0) {
                if (sum(disp.all) > 0) {
                  replace.9<-sample(disp.all, replace=FALSE, prob=dispprob,1)
                  pos.all.9<-which(replace.9==disp.all)
                  pos.9<-sample(pos.all.9,1)
                  disp.all[pos.9]<-0
                  if (pos.9 <= 9) {
                    disp.ad.vec[pos.9]<-0
                  }
                } else{
                  if (sum(disp.all) == 0) {
                    replace.9<-0
                  }   
                }
                if (replace.9 == 1) {
                  ns.9[1:10,k+1]<-0
                  ns.9[11,k+1]<-1
                  alpha.no.9 <- 1
                } else {
                  if (replace.9 > 1) {
                    ns.9[1:9,k+1]<-0
                    ns.9[10,k+1]<-replace.9-1
                    ns.9[11,k+1]<-1
                    death.9<-1
                    alpha.no.9<-1
                  } else {
                    if (replace.9 == 0) {
                      ns.9[1:11,k+1]<-0
                      death.9<-1
                      alpha.no.9<-0
                    }
                  }
                }
                last.9<-k
                
                
                
              }
            }
          }
        }
      }
    }
    
#these all just calcualte a bunch of outputs
    pck.sizes.1 <- apply(ns.1[1:11,],2,sum)
    pck.sizes.2 <- apply(ns.2[1:11,],2,sum)
    pck.sizes.3 <- apply(ns.3[1:11,],2,sum)
    pck.sizes.4 <- apply(ns.4[1:11,],2,sum)
    pck.sizes.5 <- apply(ns.5[1:11,],2,sum)
    pck.sizes.6 <- apply(ns.6[1:11,],2,sum)
    pck.sizes.7 <- apply(ns.7[1:11,],2,sum)
    pck.sizes.8 <- apply(ns.8[1:11,],2,sum)
    pck.sizes.9 <- apply(ns.9[1:11,],2,sum)
    pck.sizes.1[pck.sizes.1==0]<-NA
    pck.sizes.2[pck.sizes.2==0]<-NA
    pck.sizes.3[pck.sizes.3==0]<-NA
    pck.sizes.4[pck.sizes.4==0]<-NA
    pck.sizes.5[pck.sizes.5==0]<-NA
    pck.sizes.6[pck.sizes.6==0]<-NA
    pck.sizes.7[pck.sizes.7==0]<-NA
    pck.sizes.8[pck.sizes.8==0]<-NA
    pck.sizes.9[pck.sizes.9==0]<-NA
    pck.vec.1[j]<-mean(pck.sizes.1, na.rm=TRUE)
    pck.vec.sd.1[j]<-sd(pck.sizes.1, na.rm=TRUE)
    pups.vec.1[pups.vec.1==0]<-NA
    pup.sd.1<-sd(pups.vec.1)
    pck.vec.2[j]<-mean(pck.sizes.2, na.rm=TRUE)
    pck.vec.sd.2[j]<-sd(pck.sizes.2, na.rm=TRUE)
    pups.vec.2[pups.vec.2==0]<-NA
    pup.sd.2<-sd(pups.vec.2)
    pck.vec.3[j]<-mean(pck.sizes.3, na.rm=TRUE)
    pck.vec.sd.3[j]<-sd(pck.sizes.3, na.rm=TRUE)
    pups.vec.3[pups.vec.3==0]<-NA
    pup.sd.3<-sd(pups.vec.3)
    pck.vec.4[j]<-mean(pck.sizes.4, na.rm=TRUE)
    pck.vec.sd.4[j]<-sd(pck.sizes.4, na.rm=TRUE)
    pups.vec.4[pups.vec.4==0]<-NA
    pup.sd.4<-sd(pups.vec.4)
    pck.vec.5[j]<-mean(pck.sizes.5, na.rm=TRUE)
    pck.vec.sd.5[j]<-sd(pck.sizes.5, na.rm=TRUE)
    pups.vec.5[pups.vec.5==0]<-NA
    pup.sd.5<-sd(pups.vec.5)
    pck.vec.6[j]<-mean(pck.sizes.6, na.rm=TRUE)
    pck.vec.sd.6[j]<-sd(pck.sizes.6, na.rm=TRUE)
    pups.vec.6[pups.vec.6==0]<-NA
    pup.sd.6<-sd(pups.vec.6)
    pck.vec.7[j]<-mean(pck.sizes.7, na.rm=TRUE)
    pck.vec.sd.7[j]<-sd(pck.sizes.7, na.rm=TRUE)
    pups.vec.7[pups.vec.7==0]<-NA
    pup.sd.7<-sd(pups.vec.7)
    pck.vec.8[j]<-mean(pck.sizes.8, na.rm=TRUE)
    pck.vec.sd.8[j]<-sd(pck.sizes.8, na.rm=TRUE)
    pups.vec.8[pups.vec.8==0]<-NA
    pup.sd.8<-sd(pups.vec.8)
    pck.vec.9[j]<-mean(pck.sizes.9, na.rm=TRUE)
    pck.vec.sd.9[j]<-sd(pck.sizes.9, na.rm=TRUE)
    pups.vec.9[pups.vec.9==0]<-NA
    pup.sd.9<-sd(pups.vec.9)
    mean.pck[j]<-mean(c(pck.vec.9[j], pck.vec.8[j],pck.vec.7[j],pck.vec.6[j],pck.vec.5[j],pck.vec.4[j],pck.vec.3[j],pck.vec.2[j], pck.vec.1[j]))
    sd.pack[j]<-mean(c(pck.vec.sd.9[j], pck.vec.sd.8[j],pck.vec.sd.7[j],pck.vec.sd.6[j],pck.vec.sd.5[j],pck.vec.sd.4[j],pck.vec.sd.3[j],pck.vec.sd.2[j], pck.vec.sd.1[j]))
    litter.vec.1[j]<-mean(pups.vec.1, na.rm=TRUE)
    nextbreed.vec.1[nextbreed.vec.1==0]<-NA
    nbreed.vec.1[j]<-mean(nextbreed.vec.1,na.rm=TRUE)
    nbreed.vec.sd.1[j]<-sd(nextbreed.vec.1,na.rm=TRUE)
    litter.vec.2[j]<-mean(pups.vec.2, na.rm=TRUE)
    nextbreed.vec.2[nextbreed.vec.2==0]<-NA
    nbreed.vec.2[j]<-mean(nextbreed.vec.2,na.rm=TRUE)
    nbreed.vec.sd.2[j]<-sd(nextbreed.vec.2,na.rm=TRUE)
    litter.vec.3[j]<-mean(pups.vec.3, na.rm=TRUE)
    nextbreed.vec.3[nextbreed.vec.3==0]<-NA
    nbreed.vec.3[j]<-mean(nextbreed.vec.3,na.rm=TRUE)
    nbreed.vec.sd.3[j]<-sd(nextbreed.vec.3,na.rm=TRUE)
    litter.vec.4[j]<-mean(pups.vec.4, na.rm=TRUE)
    nextbreed.vec.4<-nextbreed.vec.4[nextbreed.vec.4==0]<-NA
    nbreed.vec.4[j]<-mean(nextbreed.vec.4,na.rm=TRUE)
    nbreed.vec.sd.4[j]<-sd(nextbreed.vec.4,na.rm=TRUE)
    litter.vec.5[j]<-mean(pups.vec.5, na.rm=TRUE)
    nextbreed.vec.5[nextbreed.vec.5==0]<-NA
    nbreed.vec.5[j]<-mean(nextbreed.vec.5,na.rm=TRUE)
    nbreed.vec.sd.5[j]<-sd(nextbreed.vec.5,na.rm=TRUE)
    litter.vec.6[j]<-mean(pups.vec.6, na.rm=TRUE)
    nextbreed.vec.6[nextbreed.vec.6==0]<-NA
    nbreed.vec.6[j]<-mean(nextbreed.vec.6,na.rm=TRUE)
    nbreed.vec.sd.6[j]<-sd(nextbreed.vec.6,na.rm=TRUE)
    litter.vec.7[j]<-mean(pups.vec.7, na.rm=TRUE)
    nextbreed.vec.7[nextbreed.vec.7==0]<-NA
    nbreed.vec.7[j]<-mean(nextbreed.vec.7,na.rm=TRUE)
    nbreed.vec.sd.7[j]<-sd(nextbreed.vec.7,na.rm=TRUE)
    litter.vec.8[j]<-mean(pups.vec.8, na.rm=TRUE)
    nextbreed.vec.8[nextbreed.vec.8==0]<-NA
    nbreed.vec.8[j]<-mean(nextbreed.vec.8,na.rm=TRUE)
    nbreed.vec.sd.8[j]<-sd(nextbreed.vec.8,na.rm=TRUE)
    litter.vec.9[j]<-mean(pups.vec.9, na.rm=TRUE)
    nextbreed.vec.9[nextbreed.vec.9==0]<-NA
    nbreed.vec.9[j]<-mean(nextbreed.vec.9,na.rm=TRUE)
    nbreed.vec.sd.9[j]<-sd(nextbreed.vec.9,na.rm=TRUE)
    mean.pups[j]<-mean(c(pups.vec.9,pups.vec.8,pups.vec.7,pups.vec.6,pups.vec.5,pups.vec.4,pups.vec.3,pups.vec.2,pups.vec.1))
    sd.pups[j]<-mean(c(pup.sd.9,pup.sd.8,pup.sd.7,pup.sd.6,pup.sd.5,pup.sd.4,pup.sd.3,pup.sd.2,pup.sd.1))
    mean.IBI[j]<-mean(c(nbreed.vec.9[j],nbreed.vec.8[j],nbreed.vec.7[j],nbreed.vec.6[j],nbreed.vec.5[j],nbreed.vec.4[j],nbreed.vec.3[j],nbreed.vec.2[j],nbreed.vec.1[j]),na.rm=TRUE)
    sd.IBI[j]<-mean(c(nbreed.vec.sd.9[j],nbreed.vec.sd.8[j],nbreed.vec.sd.7[j],nbreed.vec.sd.6[j],nbreed.vec.sd.5[j],nbreed.vec.sd.4[j],nbreed.vec.sd.3[j],nbreed.vec.sd.2[j],nbreed.vec.sd.1[j]),na.rm=TRUE)
    blank.vec.1[j]<-blank.1-1
    blank.vec.2[j]<-blank.2-1
    blank.vec.3[j]<-blank.3-1
    blank.vec.4[j]<-blank.4-1
    blank.vec.5[j]<-blank.5-1
    blank.vec.6[j]<-blank.6-1
    blank.vec.7[j]<-blank.7-1
    blank.vec.8[j]<-blank.8-1
    blank.vec.9[j]<-blank.9-1
    mean.blank[j]<-mean(c(blank.vec.1[j],blank.vec.2[j],blank.vec.3[j],blank.vec.4[j],blank.vec.5[j],blank.vec.6[j],blank.vec.7[j],blank.vec.8[j],blank.vec.9[j],blank.vec.1[j]))
    sd.blank[j]<-sd(c(blank.vec.1[j],blank.vec.2[j],blank.vec.3[j],blank.vec.4[j],blank.vec.5[j],blank.vec.6[j],blank.vec.7[j],blank.vec.8[j],blank.vec.9[j],blank.vec.1[j]))
    end.all.vec[j]<-min(end.vec[end.vec>1],na.rm=TRUE)
    count.packs<-count.packs[2:n.steps]
    count.all.packs[j]<-mean(count.packs[count.packs>0])
    dfall[j,]<-all.ad.packs
  }
  
  
  
#})
  
  
    
mean(mean.pck)
mean(mean.pups)
mean(sd.pups)
mean(sd.pack)
mean(sd.IBI)
mean(mean.IBI)
mean(mean.blank)
sd(mean.blank)
mean(end.all.vec, na.rm=TRUE)
sd(end.all.vec, na.rm=TRUE)
mean(count.all.packs)
sd(count.all.packs)


proc.time() - ptm





#for(i in 1:n.steps){
#  dfall[repeats+1,i]<-mean(dfall[1:repeats,i])
#}
#plot(dfall[6,])

#View(dfall[,500:600])

#View(dispdf)


#View(ns.1)

#mean(litter.vec)
#mean(nbreed.vec,na.rm=TRUE)

traceback()

View(ns.3)
