
rm(list=ls())

Rcpp::sourceCpp('AMS-cal2.cpp')
get.oc <- function(target, pE.true,pT.true,  ncohort, cohortsize, startdose=1, cutoff.eli=0.95, 
                   lfT, lfE, lfC, pEE.sample, ntrial=10)
{ 
  set.seed(6);
  ndose=length(pE.true)	
  npts = ncohort*cohortsize;
  YT=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
  YE=matrix(rep(0, ndose*ntrial), ncol=ndose); # store toxicity outcome
  N=matrix(rep(0, ndose*ntrial), ncol=ndose); # store the number of patients
  dselect = rep(0, ntrial); # store the selected dose level
  sel=rep(0,ndose);
  pts=rep(0,ndose);
  dlt=rep(0,ndose);
  eff=rep(0,ndose);
  tox=0;
  acr=0;
  exc=0;
  
  ################## simulate trials ###################
  for(trial in 1:ntrial)
  {
    ptm <- proc.time()
    yT<-yE<-rep(0, ndose);    ## the number of DLT at each dose level
    n<-rep(0, ndose);    ## the number of patients treated at each dose level
    earlystop=0;         ## indiate whether the trial terminates early
    d=startdose;         ## starting dose level
    elimi = rep(0, ndose);  ## indicate whether doses are eliminated
    likeliT<-matrix(rep(1,NN*(ndose+1)),ncol=ndose+1)
    likeliE<-matrix(rep(1,NN*ndose),ncol=ndose)
    likeliC<-matrix(rep(1,NN*ndose),ncol=ndose)
    likT<-rep(1,ndose+1)
    likE<-rep(1,ndose)
    UT<-rep(0,ndose)
    UE<-rep(0,ndose) 
    for(i in 1:ncohort)  
    { 
      ### generate toxicity outcome
      wT = sum(runif(cohortsize)<pT.true[d])
      yT[d] = yT[d] + wT;
      wE = sum(runif(cohortsize)<pE.true[d])
      yE[d] = yE[d] + wE;
      n[d] = n[d] + cohortsize;
      nc = n[d]/cohortsize;
      temp = post_cal2(likeliT,likeliE,likeliC,lfT[,(d-1)*(cohortsize+1)+wT+1,],
                      lfE[,(d-1)*(cohortsize+1)+wE+1,],lfC[,(d-1)*(cohortsize+1)+wE+1,],ndose);
      likeliT = temp$likeliT;
      likeliE = temp$likeliE;
      likeliC = temp$likeliC;
      posT = temp$posT;
      if(n[d]>=12 & max(yE[d]/(n[d]))>=0.25){
        posE = temp$posC;
      } else {posE = temp$posE;}
      # likeliT<-likeliT*lfT[,(d-1)*(cohortsize+1)+wT+1,]
      # likeliE<-likeliE*lfE[,(d-1)*(cohortsize+1)+wE+1,]
      # likT<-apply(likeliT,2,mean)
      # likE<-apply(likeliE,2,mean)
      # posT<-likT/sum(likT)
      # posE<-likE/sum(likE)
      cumposT<-cumsum(posT)
      dT_opt<-max(min(which(cumposT>0.85))-1,1)
      #dT_opt<-max(which.max(posT)-1,1)
      dE_opt<-which.max(posE)
      d_opt<-min(dE_opt,dT_opt)
      
      if ((1-pbeta(target,1+yT[d],n[d]-yT[d]+1))>0.95) {elimi[d:ndose]<-1
      if (elimi[1]==1) {earlystop=1; break;}}
      if (pbeta(0.25,1+yE[d],1+n[d]-yE[d])>0.9){elimi[d]<-1}
      if (sum(elimi==1)==ndose){earlystop=1;break;}
      
      if (sum(yT+yE)==0 & d<ndose){
        d<-d+1
      } else{
        
        if (d_opt>d & sum(elimi[(d+1):ndose]==0)>0){
          d<-d+which(elimi[(d+1):ndose]==0)[1]
        } else if (d_opt<d & sum(elimi[1:(d-1)]==0)>0){
          d<-max(which(elimi[1:(d-1)]==0))
        }
        
      }            
      if (sum(elimi==1)==ndose){earlystop=1;break;}
      if(elimi[d]==1){earlystop=1;break}           
      
    }
    
    YT[trial,]=yT;
    YE[trial,]=yE;
    N[trial,]=n;
    if (earlystop==0){
      posE<-posteriorHE(yE,n,pEE.sample) 
      dE_opt<-which.max(posE)
      dE_opt<-min(which(posE>=(max(posE)-0.1)))
      dE_opt<-min(which(posE>=(max(posE)-0.05)))		
      d_opt<-min(dE_opt,dT_opt)  		
      dselect[trial]=d_opt
      sel[dselect[trial]]=sel[dselect[trial]]+1/ntrial*100
    } else {dselect[trial]<-99}	
    pts<-pts+n/ntrial
    dlt<-dlt+yT/ntrial
    eff<-eff+yE/ntrial
    #print(c(trial,proc.time()-ptm))
  }	
  pts<-round(pts,1)
  dlt<-round(dlt,1)
  #print(sel)
  #print(pts) 
  results=list(sel=sel,pts=pts,dlt=dlt,eff=eff)
  return(results)	
}

sim.bams = function(act.pp, tox.pp, p0, n, M, nr) {
  safe.doses = which(tox.pp<=p0)
  tru.mtd = ifelse(length(safe.doses)>0, max(safe.doses), 0)
  max.act = max(act.pp)
  tru.mfad = min(which(act.pp==max.act))
  if (length(safe.doses)==0) tru.obd = 0 else {
    act.pp.safe = act.pp[safe.doses]
    max.act.safe = max(act.pp.safe)
    tru.obd = min(which(act.pp.safe==max.act.safe))
  }
  d = length(act.pp)
  sim.rst = get.oc(target=p0, pE.true=act.pp, pT.true=tox.pp, ncohort=M, cohortsize=n, startdose=1,
                  lfT=lfT, lfE=lfE, lfC=lfC, pEE.sample=pEE.sample, cutoff.eli=0.95, ntrial=nr)
  pct.obd = ifelse(tru.obd==0, 0, sim.rst$sel[tru.obd])
  N.ave = sum(sim.rst$pts)
  prop.trt = sim.rst$pts/N.ave
  pct.mtd.over = 100*ifelse(tru.mtd==d, 0, sum(prop.trt[(tru.mtd+1):d]))
  pct.dlt = 100*sum(prop.trt*tox.pp)
  pct.act = 100*sum(prop.trt*act.pp)
  c(pct.obd, pct.mtd.over, pct.dlt, pct.act, sim.rst$sel, N.ave, nr)
}



target = 0.25 # 0.28 # 
cohortsize = 3
ncohort = 12

ndose = 6;
NN<-50000
pT.sample<-matrix(0,NN,ndose*(ndose+1));

pT.sample[,1]<-runif(NN,target,1)
for (i in 2:ndose){
  pT.sample[,i]<-runif(NN,pT.sample[,i-1],1)
}

for (i in 1:ndose){
  pT.sample[,i+i*ndose]<-runif(NN,0,target)
  if (i>1){
    for (j in (i-1):1){
      pT.sample[,j+i*ndose]<-runif(NN,0,pT.sample[,j+1+i*ndose])
    }
  }
  if (i<ndose){
    for (j in (i+1):ndose){
      pT.sample[,j+i*ndose]<-runif(NN,sapply(pT.sample[,j-1+i*ndose],function(x) max(x,target)),1)
    }
  }
  
}


pE.sample<-matrix(0,NN,ndose*ndose);

for (i in 1:ndose){
  pE.sample[,i+(i-1)*ndose]<-runif(NN,0.25,1)
  if (i>1){
    for (j in (i-1):1){
      pE.sample[,j+(i-1)*ndose]<-runif(NN,0,pE.sample[,j+1+(i-1)*ndose])
    }
    #pE.sample[,1:i+(i-1)*ndose]<-t(apply(pE.sample[,1:i+(i-1)*ndose],1,sort))
  }
  if (i<ndose){
    for (j in (i+1):ndose){
      pE.sample[,j+(i-1)*ndose]<-runif(NN,0,pE.sample[,j-1+(i-1)*ndose])
    }
    #pE.sample[,i:ndose+(i-1)*ndose]<-t(apply(pE.sample[,i:ndose+(i-1)*ndose],1,sort,decreasing=T))
  }
  
}


pC.sample<-matrix(0,NN,ndose*ndose);

for (i in 1:ndose){
  pC.sample[,i+(i-1)*ndose]<-runif(NN,0.5,1)
  if (i>1){
    for (j in (i-1):1){
      pC.sample[,j+(i-1)*ndose]<-runif(NN,0,pC.sample[,j+1+(i-1)*ndose])
    }
    #pE.sample[,1:i+(i-1)*ndose]<-t(apply(pE.sample[,1:i+(i-1)*ndose],1,sort))
  }
  if (i<ndose){
    for (j in (i+1):ndose){
      pC.sample[,j+(i-1)*ndose]<-runif(NN,0,pC.sample[,j-1+(i-1)*ndose])
    }
    #pE.sample[,i:ndose+(i-1)*ndose]<-t(apply(pE.sample[,i:ndose+(i-1)*ndose],1,sort,decreasing=T))
  }
  
}

pEE.sample<-matrix(0,NN,ndose*ndose);

for (i in 1:ndose){
  pEE.sample[,i+(i-1)*ndose]<-runif(NN,0,1)
  if (i>1){
	for(j in (i-1):1){
	  temp<-rbinom(NN,1,0.7)
      pEE.sample[,j+(i-1)*ndose]<-temp*runif(NN,0,pEE.sample[,j+1+(i-1)*ndose])+(1-temp)*pEE.sample[,j+1+(i-1)*ndose]
	  }
  }
  if (i<ndose){
    for (j in (i+1):ndose){
	  temp<-rbinom(NN,1,0.7)
      pEE.sample[,j+(i-1)*ndose]<-temp*runif(NN,0,pEE.sample[,j-1+(i-1)*ndose])+(1-temp)*pEE.sample[,j-1+(i-1)*ndose]
    }
    #pEE.sample[,i:ndose+(i-1)*ndose]<-t(apply(pEE.sample[,i:ndose+(i-1)*ndose],1,sort,decreasing=T))
  }
  
}

lfT<-array(dim=c(NN,ndose*(cohortsize+1),ndose+1))

for (i in 1:(ndose+1)){
  for (j in 1:ndose){
    for (k in 0:cohortsize){
      lfT[,(j-1)*(cohortsize+1)+k+1,i]<-pT.sample[,j+(i-1)*ndose]^k*(1-pT.sample[,j+(i-1)*ndose])^(cohortsize-k)
    }
  }
}

lfE<-array(dim=c(NN,ndose*(cohortsize+1),ndose))

for (i in 1:ndose){
  for (j in 1:ndose){
    for (k in 0:cohortsize){
      lfE[,(j-1)*(cohortsize+1)+k+1,i]<-pE.sample[,j+(i-1)*ndose]^k*(1-pE.sample[,j+(i-1)*ndose])^(cohortsize-k)
    }
  }
}

lfC<-array(dim=c(NN,ndose*(cohortsize+1),ndose))

for (i in 1:ndose){
  for (j in 1:ndose){
    for (k in 0:cohortsize){
      lfC[,(j-1)*(cohortsize+1)+k+1,i]<-pC.sample[,j+(i-1)*ndose]^k*(1-pC.sample[,j+(i-1)*ndose])^(cohortsize-k)
    }
  }
}

p0 = target # 0.25 # 0.28 # 
tox.pp = c(0.01, 0.02, 0.04, 0.10, 0.20, 0.40) # c(0.15, 0.30, 0.45, 0.60, 0.75, 0.90) # c(0.05, 0.15, 0.30, 0.45, 0.60, 0.75) # 
d = length(tox.pp)
act.pp.mx = rbind(c(0.10, 0.20, 0.30, 0.33, 0.35, 0.36),
                  c(0.10, 0.20, 0.30, 0.30, 0.30, 0.35),
                  c(0.10, 0.20, 0.25, 0.30, 0.30, 0.35),
                  c(0.10, 0.20, 0.25, 0.30, 0.30, 0.30),
                  c(0.10, 0.20, 0.30, 0.30, 0.30, 0.30),
                  c(0.10, 0.20, 0.30, 0.30, 0.28, 0.25),
                  c(0.10, 0.20, 0.30, 0.28, 0.25, 0.20),
                  c(0.10, 0.20, 0.30, 0.28, 0.30, 0.32),
                  c(0.10, 0.20, 0.30, 0.25, 0.35, 0.35))
snrs = 1:nrow(act.pp.mx) # activity scenarios
ns = length(snrs)
M = 2*d # d # 
n = cohortsize
nr = 10^3
sim.smry = matrix(NA, ns, d+6)
file.out = "sim rst 230710a BAMS M=12 med-tox.csv"

for (s in 1:ns) {
  print(paste(s, date()))
  snr = snrs[s]
  act.pp = act.pp.mx[snr,]
  sim.smry[s,] = sim.bams(act.pp, tox.pp, p0, n, M, nr)
}

print(sim.smry)
write.csv(sim.smry, file=file.out)

# end of program