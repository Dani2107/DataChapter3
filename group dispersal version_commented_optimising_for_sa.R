#Dummy individual based model
library(pracma)
library(msm)
library(matrixStats)

#To do:
# add burn in time - set temp to average for first 10 years or whatever
# run it for 10 years with average temperature then take pack numbers and dispersers from that then run from those
# have dispersal probability as an output
# track how long packs last?

# demographic functions
#function that calculates survival probablility of adults
ad_surv<-function(temp,size,a,b,c){
  u<-a*exp(b*temp+c*size)
  u<-1-u
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
  u<-a*exp(b*size)
  return(u)}

#function that calculates number of pups
pup_no<-function(size,a,b){
  u<-exp(a+b*size)
  u<-ceiling(u)
  return(u)
}

#function for IBI
IBI<-function(temp,pups,a,b,c){
  u<-a+b*temp+c*pups
  u<-round(u)
  return(u)
}

#functions in loop
#alpha death - calculates if the alpha dies or not
alpha.die<-function(a,b){
  c<-runif(a,0,1) < b
}

#sum adult death - calculates the total number of adults (excluding the alpha) that die
all_death<-function(a,b){
  length_a = length(a)
  result = numeric(length_a)
  for(count in 1:length_a){
     c<-(runif(a[count],0,1) < b[count])
     result[count] <- length(c[c==FALSE])
  }
  return(result)
}


#sum adult disp - calculates the number of adults that disperse
adult_disp<-function(a,b){
  d <- matrix(0,nrow=1,ncol=9)

  for(counter in 1:length(a)){
    c<-(runif(a[counter],0,1) < b[counter])
    if(length(c[c==TRUE])>=1){
      d[counter]<-sum(round(rtnorm((length(c[c==TRUE])),3,1,lower=0)))
    }
}
 return(d)
}

#calculate lambda
lambda_func<-function(a,b){
  a<-mean(a[b-12:b])/mean(a[12:24])
}


run_simulation <- function(ad.surv.intercept,ad.surv.temp,
                           ad.surv.pack,disp.intercept,disp.pack.size,IBI.intercept,IBI.tempcoef,IBI.pups,
                           pup.intercept,pup.pack.size,juv.surv.intercept,juv.surv.temp,juv.surv.litter,
                           meantemp,actual.meantemp,number.of.packs,n.steps,repeats,burntime)
{
  #n.steps  - number of timesteps to run the model for ech time
  #repeats  - number of times to repeat the loop
  dfall<-array(0,c(repeats+1,n.steps-(burntime-1))) #number of dogs in the population t each timestep across runs
  #vectors of pack size, mean litter size, IBI, n timesteps empty, sd of pack size, IBI and pup number
  pup.sd <- matrix(0,nrow=1,ncol=number.of.packs)
  pck.vec <- matrix(0,nrow=repeats,ncol=number.of.packs)
  pck.vec.sd <-matrix(0,nrow=repeats,ncol=number.of.packs)
  litter.vec <- matrix(0,nrow=repeats,ncol=number.of.packs)
  nbreed.vec <- matrix(0,nrow=repeats,ncol=number.of.packs)
  nbreed.vec.sd <- matrix(0,nrow=repeats,ncol=number.of.packs)
  blank.vec <-matrix(0,nrow=repeats,ncol=number.of.packs)
  alpha.no <-  matrix(1,nrow=1,ncol=number.of.packs) #number of alphas in a particular matrix (pack)
  mean.pck<-sd.pack<-mean.pups<-mean.IBI<-sd.IBI<-mean.blank<-sd.blank<-end.all.vec<-count.all.packs<-rep(0,repeats)
  dispdf<-array(0,c(number.of.packs*2,n.steps+1))
  disp.vec<-rep(0,repeats)
  lambda.vec<-rep(0,repeats)
  all.dogs.vec<-rep(0,repeats)
  pack.age.vec<-rep(0,repeats)
  sd.pack.age<-rep(0,repeats)
  sd.pups<-rep(0,repeats)
  sd.dogs<-rep(0,repeats)


  for(j in 1:repeats){
    # want to store the matrices, and the vector of population sizes
    ns <- array(0,dim=c(number.of.packs,11,n.steps+1))
    agepack <- matrix(NA,nrow=n.steps+1,ncol=number.of.packs)
    agetrack <- matrix(0,nrow=1,ncol=number.of.packs)
    # need to start with a population structure at time i=1
    ns[,10,1] <- rpois(number.of.packs,3)
    ns[,11,1] <- 1
    #sets my litter size 1st run
    meanlitter<-round(rnorm(1,4,0.75))
    #sets when breeding will occurr first round
    nextbreed <- matrix(3,nrow=number.of.packs,ncol=1)
    k<-0 #counts from 0 each run
    k.vec <- temp.vec<-rep(0,n.steps) #stores k and temp in a vector
    #vector for pup number, pack number, IBI, number of packs, when there are no dogs left and the total number of dogs
    count.packs<-end.vec<-all.ad.packs<-rep(0,n.steps)
    pups.vec <- matrix(0,nrow=n.steps,ncol=number.of.packs)
    pack.vec <- matrix(0,nrow=n.steps,ncol=number.of.packs)
    nextbreed.vec <-matrix(0,nrow=n.steps,ncol=number.of.packs)

    disp.ad.vec<-disp.ad.prev<-rep(0,number.of.packs) #vector of dispersers and vector dispersers move into after 1 timestep
    #disp.ad.prev2<-rep(0,9) in case want to stay for 3 months
    #pack dispersal, the final run if all dog die, steps where there are no dogs and number of dogs in the pack
    blank <- matrix(0,nrow=1,ncol=number.of.packs)
    packnumb <- matrix(0,nrow=1,ncol=number.of.packs)
    packdisp <- matrix(0,nrow=1,ncol=number.of.packs)
    last <- matrix(0,nrow=1,ncol=number.of.packs)
    
      for (i in 1:n.steps){
      k<-k+1
      k.vec[i]<-k
      

      packnumb = rowSums(ns[,1:11,k-1])
      if(k>=burntime){
      blank[packnumb == 0] <- blank[packnumb == 0] + 1
      agetrack[1,1:9]<-agetrack[1,1:9]+1
      agepack[i,packnumb == 0]<-agetrack[1,packnumb == 0]
      agetrack[1,packnumb == 0]<-0}
      all.of.packs<-sum(packnumb)
      all.ad.packs[i]<-all.of.packs
      if(all.of.packs==0){
        end.vec[i]<-k
      }

      count.packs[i]<-sum(packnumb>0)
      #death
      death<-ns[,11,k] #if there is an alpha death[n] = 1, otherwise 0
      #define the size of the pack
      pack.size <- rowSums(ns[,1:11,i])
      pack.vec[i,]<-pack.size
      #calculates litter size first run
      litter.size <- pups <- pup_no(pack.size,pup.intercept,pup.pack.size)
      pups.vec[i,]<-pups
      # generate temperature
      if(j>repeats/variables*13&&j<=repeats/variables*14){
      if(k<=burntime){
        temperature <- rnorm(1,actual.meantemp*1.01,1)}
      else {
        temperature <- rnorm(1,meantemp*1.01,1)
      }
      } else {
        if(k<=burntime){
          temperature <- rnorm(1,actual.meantemp,1)}
        else {
          temperature <- rnorm(1,meantemp,1)
        }
      }
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
      if(j>=0&&j<=repeats/variables){
        surv.age.a <- ad_surv(temperature,pack.size,ad.surv.intercept*1.01,ad.surv.temp,ad.surv.pack)
      } else if (j>repeats/variables&&j<=repeats/variables*2) {
        surv.age.a <- ad_surv(temperature,pack.size,ad.surv.intercept,ad.surv.temp*1.01,ad.surv.pack)
      } else if (j>repeats/variables*2&&j<=repeats/variables*3) {
        surv.age.a <- ad_surv(temperature,pack.size,ad.surv.intercept,ad.surv.temp,ad.surv.pack*1.01)
      } else {
        surv.age.a <- ad_surv(temperature,pack.size,ad.surv.intercept,ad.surv.temp,ad.surv.pack)}
      #number of non alpha adults
      ad.numb<-ns[,10,k]
      #number of adults that die
      death.adult = all_death(ad.numb,surv.age.a)
      #new number of adults after death
      ns[,10,k+1] <- (ns[,10,k]-death.adult)
      #dispersal
      if (j>repeats/variables*3&&j<=repeats/variables*4) {
        dispersal<-ad_disp(pack.size,disp.intercept*1.01,disp.pack.size)
      } else if (j>repeats/variables*4&&j<=repeats/variables*5){
        dispersal<-ad_disp(pack.size,disp.intercept,disp.pack.size*1.01)
      } else {
        dispersal<-ad_disp(pack.size,disp.intercept,disp.pack.size)}
      #updates ad.numb to post death figures
      ad.numb<-ns[,10,k+1]
      #calculates the number of dispersing adults
      disp.adult = matrix(0,nrow=1,ncol=number.of.packs)
      disp.adult <-adult_disp(ad.numb,dispersal)
      disp.adult[disp.adult > ad.numb] = ad.numb[disp.adult > ad.numb]

      #adds dispersers to the dispersal vector
      disp.ad.vec <- disp.adult + packdisp
      #object that gets filled with any pack dispersal numbers resulting from alpha death
      packdisp <- matrix(0,nrow=1,ncol=number.of.packs)
      #updates matrix with number of adults after dispersal
      ns[,10,k+1] <- (ns[,10,k+1] - disp.adult)
      #set juvenile survival
      if (j>repeats/variables*7&&j<=repeats/variables*8){
        surv.juv = juv_surv(temperature,litter.size,juv.surv.intercept*1.01,juv.surv.temp,juv.surv.litter)
      } else if (j>repeats/variables*8&&j<=repeats/variables*9) {
        surv.juv = juv_surv(temperature,litter.size,juv.surv.intercept,juv.surv.temp*1.01,juv.surv.litter)
      } else if (j>repeats/variables*9&&j<=repeats/variables*10) {
        surv.juv = juv_surv(temperature,litter.size,juv.surv.intercept,juv.surv.temp,juv.surv.litter*1.01)
      } else {
        surv.juv = juv_surv(temperature,litter.size,juv.surv.intercept,juv.surv.temp,juv.surv.litter)
      }
      #this calculates number of dogs that die each timestep.
      for(loop_counter in 1:number.of.packs){
        #gives which age class has juveniles in it
        nonzeroes.juv <- (which(ns[loop_counter,0:number.of.packs,k]>=1))
        #number of juvenile age classes with dogs in
        length.juv <- length(nonzeroes.juv)
        if(length.juv == 1) {
        nonzeroes.juv.first <- (nonzeroes.juv[1])
        if(nonzeroes.juv.first < 9) {
          nonzeroes.juv.n.first <- (ns[loop_counter,nonzeroes.juv.first,k])
          death.juv.first <- all_death(nonzeroes.juv.n.first,surv.juv[loop_counter])
          ns[loop_counter,nonzeroes.juv.first + 1,k+1]<-(ns[loop_counter,nonzeroes.juv.first,k] - death.juv.first)
        } else {
          nonzeroes.juv.n.first <- (ns[loop_counter,nonzeroes.juv.first,k])
          death.juv.first <- all_death(nonzeroes.juv.n.first,surv.juv[loop_counter])
          ns[loop_counter,nonzeroes.juv.first+1,k+1]<-((ns[loop_counter,nonzeroes.juv.first,k]-death.juv.first) + ns[loop_counter,10,k])
        }
      }
      if(length.juv > 1) {
        nonzeroes.juv.second <- (nonzeroes.juv[2])
        nonzeroes.juv.n.second <- (ns[loop_counter,nonzeroes.juv.second,k])
        death.juv.second <- all_death(nonzeroes.juv.n.second,surv.juv[loop_counter])
        ns[loop_counter,10,k+1]<-(ns[loop_counter,nonzeroes.juv.second,k]-death.juv.second)+ns[loop_counter,10,k]
        nonzeroes.juv.first <- (nonzeroes.juv[1])
        nonzeroes.juv.n.first <- (ns[loop_counter,nonzeroes.juv.first,k])
        death.juv.first <- all_death(nonzeroes.juv.n.first,surv.juv[loop_counter])
        ns[loop_counter,nonzeroes.juv.first+1,k+1]<-(ns[loop_counter,nonzeroes.juv.first,k]-death.juv.first)
        }
      }
      # now do the probability of recruitment - ie do they give birth or not, varies with temp
      pb =  matrix(0,nrow=number.of.packs,ncol=1)
      pb[i == nextbreed] = 1

      # Pick a the next step (i) when breeding occurs.
      litterdiff = matrix(0,nrow=number.of.packs,ncol=1)
      nextbreed.sum = matrix(0,nrow=number.of.packs,ncol=1)

      tempdiff <- breedtemp - actual.meantemp
      pb_true <- pb == 1
      num_pb_true = sum(pb_true)
      litterdiff[pb_true] <- pups[pb_true] - meanlitter
      nextbreed.sum[pb_true] <- nextbreed[pb_true]
      if (j>repeats/variables*10&&j<=repeats/variables*11){
        nextbreed[pb_true] <- round(rnorm(num_pb_true, mean = (IBI(breedtemp,pups,IBI.intercept*1.01,IBI.tempcoef,IBI.pups)), 1) + i)
      } else if (j>repeats/variables*11&&j<=repeats/variables*12) {
        nextbreed[pb_true] <- round(rnorm(num_pb_true, mean = (IBI(breedtemp,pups,IBI.intercept,IBI.tempcoef*1.01,IBI.pups)), 1) + i) 
      } else if (j>repeats/variables*12&&j<=repeats/variables*13) {
        nextbreed[pb_true] <- round(rnorm(num_pb_true, mean = (IBI(breedtemp,pups,IBI.intercept,IBI.tempcoef,IBI.pups*1.01)), 1) + i)
      } else {
        nextbreed[pb_true] <- round(rnorm(num_pb_true, mean = (IBI(breedtemp,pups,IBI.intercept,IBI.tempcoef,IBI.pups)), 1) + i)
      }
      nextbreed.vec[i,pb_true] <- nextbreed[pb_true] - nextbreed.sum[pb_true]
      # now number of pups produced
      if (j>repeats/variables*5&&j<=repeats/variables*6){
        pups <- pup_no(pack.size,pup.intercept*1.01,pup.pack.size)
      } else if (j>repeats/variables*6&&j<=repeats/variables*7) {
        pups <- pup_no(pack.size,pup.intercept,pup.pack.size*1.01)
      } else {
        pups <- pup_no(pack.size,pup.intercept,pup.pack.size)
      }
      litter.size <- pb * pups
      which_are_true <- (ns[,11,k] == 1)
      ns[which_are_true,1,k+1] <- pb[which_are_true]*pups[which_are_true]

      #alpha death dynamic
      # now dominant suvival
      #this sets whether the alpha has died (0) or not (1) and assigns it to 'death'
      which_die <- alpha.no == 1
      if (j>repeats/variables*14&&j<=repeats/variables*15) {
        alpha.no[which_die] <- death[which_die] <- alpha.die(sum(which_die),(surv.age.a[which_die]-alphadieprob*1.01))
      } else {
        alpha.no[which_die] <- death[which_die] <- alpha.die(sum(which_die),(surv.age.a[which_die]-alphadieprob))
      }
      ns[,11,k+1] <- alpha.no
      #this determines what will happen if the alpha has died. Either the whole pack disperses
      # or another adult takes over. If there are juveniles present during whole pack dispersal they die
      replace <- matrix(0,nrow=number.of.packs,ncol=1)
      for (loop_counter in 1:number.of.packs){
      if (death[loop_counter] == 0) {
        dth <- runif(1,0,100)
        if (dth > 40) {
          ns[loop_counter,11,k+1] <- 0
          last[loop_counter] <- k
          #if the whole pack has dispersed + the matrix is empty then this picks which dispersers
          #fill the matrix (form a new pack). One becomes an alpha. if theres >1 the others are
          #non alpha adults
          if (sum(disp.all) > 0) {
            replace[loop_counter]<-sample(disp.all, replace=FALSE, prob = dispprob,1)
            pos.all <- which(replace[loop_counter]==disp.all)
            pos<-sample(pos.all,1)
            disp.all[pos]<-0
            if (pos <= number.of.packs) {
              disp.ad.vec[pos]<-0
            }
          } else {
            if (sum(disp.all) == 0) {
              replace[loop_counter]<-0
            }
          }
          if (replace[loop_counter] == 1) {
            ns[loop_counter,1:10,k+1]<-0
            ns[loop_counter,11,k+1]<-1
            death[loop_counter]<-1
            alpha.no[loop_counter]<-1
          } else {
            if (replace[loop_counter] > 1) {
              ns[loop_counter,1:9,k+1]<-0
              ns[loop_counter,10,k+1]<-replace[loop_counter]-1
              ns[loop_counter,11,k+1]<-1
              death[loop_counter]<-1
              alpha.no[loop_counter]<-1
            } else {
              if (replace[loop_counter] == 0) {
                ns[loop_counter,1:11,k+1]<-0
                death[loop_counter]<-1
                alpha.no[loop_counter]<-0
              }
            }
          }
          packdisp[loop_counter]<-disp.adult[loop_counter]+sum(ns[loop_counter,9:10,i])

        } else {

          while (death[loop_counter] == 0) {
            if (ns[loop_counter,10, k] >= 1) {

              ns[loop_counter,10, k+1] <- ns[loop_counter,10, k] - 1
              death[loop_counter] <- 1
              alpha.no[loop_counter] <-1
              ns[loop_counter,11,k+1]<-death[loop_counter]

            } else {

              # These two statements stop the loop going on forever.
              if (ns[loop_counter,10, k] == 0) {
                if (sum(disp.all) > 0) {
                  replace[loop_counter]<-sample(disp.all, replace=FALSE, prob=dispprob,1)
                  pos.all<-which(replace[loop_counter]==disp.all)
                  pos<-sample(pos.all,1)
                  disp.all[pos]<-0
                  if (pos <= number.of.packs) {
                    disp.ad.vec[pos]<-0
                  }
                } else{
                  if (sum(disp.all) == 0) {
                    replace[loop_counter]<-0
                  }
                }
                if (replace[loop_counter] == 1) {
                  ns[loop_counter,1:10,k+1]<-0
                  ns[loop_counter,11,k+1]<-1
                  alpha.no[loop_counter] <- 1
                } else {
                  if (replace[loop_counter] > 1) {
                    ns[loop_counter,1:9,k+1]<-0
                    ns[loop_counter,10,k+1]<-replace[loop_counter]-1
                    ns[loop_counter,11,k+1]<-1
                    death[loop_counter]<-1
                    alpha.no[loop_counter]<-1
                  } else {
                    if (replace[loop_counter] == 0) {
                      ns[loop_counter,1:11,k+1]<-0
                      death[loop_counter]<-1
                      alpha.no[loop_counter]<-0
                    }
                  }
                }
                last[loop_counter]<-k

              }
            }
          }
        }
      }
      }

      }

#these all just calculate a bunch of outputs

    pck.sizes = apply(ns[,1:11,],1,colSums)
    pck.sizes[pck.sizes==0]<-NA
    

    pups.vec[pups.vec==0] <- NA
    pup.sd<-colSds(pups.vec[burntime:n.steps,])

    pck.vec[j,]<-colMeans(pck.sizes[burntime:n.steps,], na.rm=TRUE)
    pck.vec.sd[j,]<-colSds(pck.sizes[burntime:n.steps,], na.rm=TRUE)

    mean.pck[j]<-mean(pck.vec[j,])
    sd.pack[j]<-mean(pck.vec.sd[j,])

    litter.vec[j,] <- colMeans(pups.vec[burntime:n.steps,], na.rm=TRUE)

    nextbreed.vec[nextbreed.vec==0]<-NA
    nbreed.vec[j,] <- colMeans(nextbreed.vec[burntime:n.steps,],na.rm=TRUE)
    nbreed.vec.sd[j,] <- colSds(nextbreed.vec[burntime:n.steps,],na.rm=TRUE)

    mean.pups[j]<-mean(pups.vec[burntime:n.steps,])
    mean.IBI[j]<-mean(nbreed.vec[j,],na.rm=TRUE)
    sd.IBI[j]<-mean(nbreed.vec.sd[j,],na.rm=TRUE)
    sd.pups[j]<-mean(pup.sd)
    sd.dogs[j]<-mean(pck.vec.sd)
    blank.vec[j,]<-blank
    mean.blank[j]<-mean(blank.vec[j,])
    end.all.vec[j]<-min(end.vec[end.vec>1],na.rm=TRUE)
    count.packs<-count.packs[burntime:n.steps]
    count.all.packs[j]<-mean(count.packs[count.packs>0])
    lambda.vec[j]<-lambda_func(all.ad.packs,n.steps)
    all.dogs.vec[j]<-mean(all.ad.packs[burntime:n.steps])
    disp.vec[j]<-mean(colSums(dispdf[,burntime:n.steps]))
    dfall[j,]<-all.ad.packs[burntime:n.steps]
    agepack[which(agepack==1)]<-NA
    pack.age.vec[j]<-mean(colMeans(agepack,na.rm=TRUE),na.rm=TRUE)/12
    sd.pack.age[j]<-mean(colSds(agepack,na.rm=TRUE),na.rm=TRUE)
  }
  
results = data.frame(mean.pck,mean.pups,mean.IBI,mean.blank,pack.age.vec,count.all.packs,lambda.vec,
                      disp.vec,all.dogs.vec,mean.blank,end.all.vec,count.all.packs,
                      sd.pack,sd.IBI, sd.pack.age,sd.pups,sd.dogs
                     ,ad.surv.intercept,ad.surv.temp,
                     ad.surv.pack,disp.intercept,disp.pack.size,IBI.intercept,IBI.tempcoef,IBI.pups,
                     pup.intercept,pup.pack.size,juv.surv.intercept,juv.surv.temp,juv.surv.litter,
                     meantemp,actual.meantemp,number.of.packs,n.steps,repeats,burntime)   

return(results)
}

summarise_results <- function(results){
  # produces summary statistics from a results data frame
  mean.pck = mean(results$mean.pck)
  mean.pups = mean(results$mean.pups)
  sd.pack = mean(results$sd.pack)
  sd.IBI = mean(results$sd.IBI)
  mean.IBI = mean(results$mean.IBI)
  mean.blank = mean(results$mean.blank)
  count.all.packs = mean(results$count.all.packs)
  
  #Get the input arguments for this simulation.
  #There should just be one set for the entire results data frame
  input_arguments = unique(results[c("ad.surv.intercept","ad.surv.temp",
                                     "ad.surv.pack","disp.intercept","disp.pack.size","IBI.intercept","IBI.tempcoef","IBI.pups",
                                     "pup.intercept","pup.pack.size","juv.surv.intercept","juv.surv.temp","juv.surv.litter",
                                     "meantemp","actual.meantemp","number.of.packs","n.steps","repeats","burntime")])
  
  summary_data <- data.frame(mean.pck,mean.pups,sd.pack,sd.IBI,mean.IBI,mean.blank,count.all.packs,
                        input_arguments
                        )
  return(summary_data)
  
}


