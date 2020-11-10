### Supplement R code FOR GRACE ET AL. 2012. GUIDELINES FOR A GRAPH-THEORETIC . . . ECOSPHERE
### Note: an error in one of the comments is corrected in this version

# CONTENTS OF FILE
# Set working directory
# Load the data
# Inspect variables
# WinBUGS work for Bayesian Analyses
#   Initial SE model.
#   SE model after pruning ns pathways.
#   Prediction equations for initial model pruned of ns pathways.
#   Examining residual relationships for missing pathways.
#   Demonstration plot for including in the journal article.
#   Final SE model
#   Prediction equations for final se model.
#   Examining residual relationships for missing pathways in final model.
#   Scenarios in Ecosphere paper.
#   Queries in lieu of standardized coefs.
# Illustration of the solution of this model using a piecewise approach and frequentist methods.


##############################
### Set your working directory
##############################
setwd("R:/83439HC - Forest and Marsh ERM/TWD - Structural Equation Modeling/Working/DATA/AcadiaNP/SEMwork/SEM-data-May2011")
setwd("F:/Documents/ms/1.Grace/Ecosphere_Guidelines-ms/RevisionForEcosphere/Appendix_A")

#####################################
### LOAD THE DATA
#####################################
surr <- c(1,0,0,1,3,1,0,0,2,2,2,3,2,3,2,0,2,0,0,2,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0)
buff <- c(0,0,0,2,3,1,0,0,3,3,3,3,2,2,3,0,2,0,0,3,0,0,1,0,0,0,1,0,0,0,0,2,0,0,0,0,0)
hyd <- c(6,1,3,3,2,5,1,2,5,4,5,5,4,3,5,1,6,0,3,3,0,0,2,0,0,0,3,0,0,0,0,2,0,0,1,1,0)
soil <- c(0,0,0,1,1,1,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
wet <- c(365.00,139.05,156.43,288.96,312.86,141.94,278.10,330.25,295.48,365.00,191.19,121.67,243.33,330.25,121.67,139.05,260.71,178.15,365.00,69.52,0.00,30.42,36.50,0.00,0.00,73.00,0.00,0.00,44.71,35.67,142.34,69.92,81.00,39.25,81.17,0.00,52.71)
cond <- c(2.113,1.899,1.730,2.127,2.583,2.038,1.816,1.995,2.392,2.510,2.653,2.737,2.692,2.278,2.718,2.117,1.817,1.820,1.577,2.115,1.946,1.764,1.725,1.624,1.705,1.813,2.363,1.646,1.920,1.893,1.415,1.990,1.781,1.943,1.726,1.641,1.858)
natr <- c(11,41,37,32,15,45,31,28,24,19,14,24,25,14,19,41,25,34,15,27,46,52,43,49,49,38,61,56,51,62,28,43,37,44,48,59,62)
sphag <- c(2,79,85,2,0,51,46,48,59,74,16,11,23,29,32,66,54,93,37,59,82,72,82,94,92,91,75,94,82,79,64,42,78,51,89,85,77)
typha <- c(0.0000,0.1239,0.0271,0.6190,1.7818,0.5346,1.0212,0.0000,0.8955,1.2672,1.3259,1.0700,1.1056,0.3044,0.4150,0.0792,0.0000,0.0000,0.0000,0.2430,0.0858,0.4001,0.0000,0.0000,0.0000,0.0000,0.7819,0.0573,0.0013,0.3884,0.0000,0.6520,0.0000,0.0000,0.0155,0.1060,0.0025)
sem.dat <- data.frame(surr,buff,hyd,soil,wet,cond,natr,sphag,typha) # sem.dat is data set used for the analyses in this paper

# variable numeric codes (used below for parameters)
# surr=1; buff=2; hyd=3; soil=4; wet=5; cond=6; natr=7; sphag=8; typha=9

####################
# INSPECT VARIABLES
####################
par(mfrow=c(3,3)); hist(surr); hist(buff); hist(hyd); hist(soil);hist(wet); hist(cond); hist(natr); hist(sphag); hist(typha)


###########################################################################
# WINBUGS WORK FOR BAYESIAN ANALYSES
###########################################################################
library(R2WinBUGS)

### SET YOUR WORKING DIRECTORY FOR STORING BAYESIAN FILES 
setwd("C:/Data/WinbugWork/GuidelinesMS") 

#############################################################################
# INITIAL SE MODEL
#############################################################################
# Creating winbugs program and saving as external file.
sink("pathsimple.init.txt")		
cat("
MODEL LR1 {
	# Likelihoods
	for(i in 1:N) {
	  buff[i] ~ dbin(buffp.hat[i],n.buff)				     # buffp is prop. of max intrusion n = max value 3
		   logit(buffp.hat[i]) <- b2.0 + b2.1*surr[i]
    	  hyd[i] ~ dbin(hydp.hat[i],n.hyd)                                   # hyd is prop. of max of transformed value; max = 6
		   logit(hydp.hat[i]) <- b3.0 + b3.1*surr[i] + b3.2*buff[i]
	  soil[i] ~ dbern(soilp.hat[i])					     # soil disturbance response - logit
		   logit(soilp.hat[i]) <- b4.0 + b4.1*buff[i]
	  wet[i] ~ dbin(wetp.hat[i],n.wet)				     # flooded response - proportional odds model
		   logit(wetp.hat[i]) <- b5.0 +b5.1*hyd[i] +b5.2*surr[i]
	  cond[i] ~ dnorm(cond.hat[i],tau.cond)			             # logged (base 1) conductivity response treated as normal
		   cond.hat[i] <- b6.0 + b6.1*soil[i] + b6.2*surr[i] + b6.3*buff[i]
	  natr[i] ~ dpois(natr.hat[i])				             # native richness response - poisson
		   log(natr.hat[i]) <- b7.0 +b7.1*wet[i] + b7.99*typha[i]
	  sphag[i] ~ dbin(sphagp.hat[i],n.sphag)			     # Lsphagnum response - proportional odds response w/ no int
		   logit(sphagp.hat[i]) <- b8.3*cond[i] + b8.6*wet[i] + b8.99*typha[i] 
	  typha[i] ~ dnorm(typha.hat[i],tau.typha)                           # Ltypha response - log normal with changepoint of psi
		typha.hat[i] <- b9.1*cond[i] + b9.2*step(cond[i]-psi.typha)*(cond[i]-psi.typha) + b9.3*wet[i]
	}
	# Priors
	b2.0 ~ dnorm(0,0.00001); b2.1 ~ dnorm(0,0.00001) 
	b3.0 ~ dnorm(0,0.00001); b3.1 ~ dnorm(0,0.00001); b3.2 ~ dnorm(0,0.00001);
	b4.0 ~ dnorm(0,0.00001); b4.1 ~ dnorm(0,0.00001);
	b5.0 ~ dnorm(0,0.00001); b5.1 ~ dnorm(0,0.00001); b5.2 ~ dnorm(0,0.00001);
	b6.0 ~ dnorm(0,0.00001); b6.1 ~ dnorm(0,0.00001); b6.2 ~ dnorm(0,0.00001); b6.3 ~ dnorm(0,0.00001);
  	tau.cond ~ dgamma(0.001,0.001)
	b7.0 ~ dnorm(0,0.00001); b7.1 ~ dnorm(0,0.00001); b7.99 ~ dnorm(0,0.00001)
	b8.3 ~ dnorm(0,0.00001); b8.6 ~ dnorm(0,0.00001); b8.99 ~ dnorm(0,0.00001)
	b9.1 ~ dnorm(0,0.00001); b9.2 ~ dnorm(0,0.00001); b9.3 ~ dnorm(0,0.00001);
	tau.typha ~ dgamma(0.001,0.001); psi.typha ~ dunif(1.5,2.5);

	# Derived Quantities
	prob.3.2.GT0 <- step(b3.2);      # these are used to calculate exact probabilities for parameters
	prob.6.3.GT0 <- step(b6.3);
	prob.7.99.GT0 <- step(b7.99);
	prob.8.99.GT0 <- step(b8.99);
	prob.9.3.GT0 <- step(b9.3);
	}
",fill=TRUE)
sink()
# end of winbugs code creation

# Creating objects to hand to winbugs.
N=37; n.buff=3; n.wet=365; n.sphag=100; n.hyd=6; #some input variables
data = list("N","n.buff","n.wet","n.sphag","n.hyd",
           "surr","buff","soil","hyd","wet","cond","natr","sphag","typha")   
parameters <- c("b2.0","b2.1","b3.0","b3.1","b3.2","b4.0","b4.1",                          
                "b5.0","b5.1","b5.2",
                "b6.0","b6.1","b6.2","b6.3",
                "b7.0","b7.1","b7.99",
                "b8.3","b8.6","b8.99",
                "b9.1","b9.2","b9.3","psi.typha",
                "prob.3.2.GT0","prob.6.3.GT0","prob.7.99.GT0","prob.8.99.GT0","prob.9.3.GT0")

inits <- function(){list(b2.0=0,b2.1=0,b3.0=0,b3.1=0,b3.2=0,b4.0=0,b4.1=0,                                
                 b5.0=0,b5.1=0,b5.2=0,
                 b6.0=0,b6.1=0,b6.2=0,b6.3=0,tau.cond=100,
                 b7.0=0,b7.1=0,b7.99=0,
                 b8.3=0,b8.6=0,b8.99=0,
		             b9.1=0,b9.2=1,b9.3=1,
                 tau.typha=100,psi.typha=2)}

# MCMC settings to hand to winbugs
 nc <- 3      # number of chains
 ni <- 3000  # number of samples for each chain      # use small number of samples to check program before final run with large number (e.g., 30,000)
 nb <- 1000   # number of samples to discard for burn in
 nt <- 1      # thinning rate

# Hand off to winbugs using "bugs" command. Note: because debug=TRUE, you will need to close winbugs before returning to R.
pathsimple.init.out <- bugs(data=data, inits=inits, parameters.to.save=parameters, model.file="pathsimple.init.txt", 
 n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE, working.directory=getwd())

# Print some basic results
print(pathsimple.init.out,digits=4)  # print parameters

# Plot posterior distributions of interest.
hist(pathsimple.init.out$sims.list$b6.1, col="blue")                  # illustration of histogram of parameter ests


#############################################################################
# SE MODEL AFTER PRUNING NS PATHWAYS
#############################################################################
# Creating winbugs program and saving as external file.
sink("pathsimple.init2.txt")		
cat("
MODEL LR1 {
	# Likelihoods
	for(i in 1:N) {
	  buff[i] ~ dbin(buffp.hat[i],n.buff)				
		logit(buffp.hat[i]) <- b2.0 + b2.1*surr[i]
    	  hyd[i] ~ dbin(hydp.hat[i],n.hyd)                      
		logit(hydp.hat[i]) <- b3.0 + b3.1*surr[i]
	  soil[i] ~ dbern(soilp.hat[i])					
		logit(soilp.hat[i]) <- b4.0 + b4.1*buff[i]
	  wet[i] ~ dbin(wetp.hat[i],n.wet)				
		logit(wetp.hat[i]) <- b5.0 +b5.1*hyd[i] +b5.2*surr[i]
	  cond[i] ~ dnorm(cond.hat[i],tau.cond)			
		cond.hat[i] <- b6.0 + b6.1*soil[i] + b6.2*surr[i] 
	  natr[i] ~ dpois(natr.hat[i])				
		log(natr.hat[i]) <- b7.0 +b7.1*wet[i] 
	  sphag[i] ~ dbin(sphagp.hat[i],n.sphag)			
		logit(sphagp.hat[i]) <- b8.3*cond[i] + b8.6*wet[i] 
	  typha[i] ~ dnorm(typha.hat[i],tau.typha)              
		typha.hat[i] <- b9.1*cond[i] + b9.2*step(cond[i]-psi.typha)*(cond[i]-psi.typha)
	}
	# Priors
	b2.0 ~ dnorm(0,0.00001); b2.1 ~ dnorm(0,0.00001) 
	b3.0 ~ dnorm(0,0.00001); b3.1 ~ dnorm(0,0.00001); 
	b4.0 ~ dnorm(0,0.00001); b4.1 ~ dnorm(0,0.00001);
	b5.0 ~ dnorm(0,0.00001); b5.1 ~ dnorm(0,0.00001); b5.2 ~ dnorm(0,0.00001);
	b6.0 ~ dnorm(0,0.00001); b6.1 ~ dnorm(0,0.00001); b6.2 ~ dnorm(0,0.00001); 
  	tau.cond ~ dgamma(0.001,0.001)
	b7.0 ~ dnorm(0,0.00001); b7.1 ~ dnorm(0,0.00001); 
	b8.3 ~ dnorm(0,0.00001); b8.6 ~ dnorm(0,0.00001);
	b9.1 ~ dnorm(0,0.00001); b9.2 ~ dnorm(0,0.00001); 
	tau.typha ~ dgamma(0.001,0.001); psi.typha ~ dunif(1.5,2.5);

	# Derived Quantities
	prob.6.1.GT0 <- step(b6.1);
	}
",fill=TRUE)
sink()
# end of winbugs code creation

# Creating objects to hand to winbugs.
N=37; n.buff=3; n.wet=365; n.sphag=100; n.hyd=6; 
data = list("N","n.buff","n.wet","n.sphag","n.hyd",
           "surr","buff","soil","hyd","wet","cond","natr","sphag","typha")   
parameters <- c("b2.0","b2.1","b3.0","b3.1","b4.0","b4.1",                          
                "b5.0","b5.1","b5.2",
                "b6.0","b6.1","b6.2",
                "b7.0","b7.1",
                "b8.3","b8.6",
                "b9.1","b9.2","psi.typha",
                "prob.6.1.GT0")

inits <- function(){list(b2.0=0,b2.1=0,b3.0=0,b3.1=0,b4.0=0,b4.1=0,                                
                 b5.0=0,b5.1=0,b5.2=0,
                 b6.0=0,b6.1=0,b6.2=0,tau.cond=100,
                 b7.0=0,b7.1=0,
                 b8.3=0,b8.6=0,
		             b9.1=0,b9.2=1,
                 tau.typha=100,psi.typha=2)}

# MCMC settings to hand to winbugs
 nc <- 3      # number of chains
 ni <- 30000  # number of samples for each chain
 nb <- 1000   # number of samples to discard for burn in
 nt <- 1      # thinning rate

# Hand off to winbugs using "bugs" command. Note: because debug=TRUE, you will need to close winbugs before returning to R.
pathsimple.init2.out <- bugs(data=data, inits=inits, parameters.to.save=parameters, model.file="pathsimple.init2.txt", 
 n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE, working.directory=getwd())

# Print some basic results
print(pathsimple.init2.out,digits=4)  # print parameters

# Plot posterior distributions of interest.
hist(pathsimple.init2.out$sims.list$b6.1, col="blue")                  # illustration of histogram of parameter ests


###########################################################
# PREDICTION EQNS FOR INITIAL SE MODEL PRUNED OF NS PATHS
###########################################################

### parameter estimates from winbugs analysis
b2.0=           -3.1012 
b2.1=            2.4392 
b3.0=           -1.5972 
b3.1=            1.1392 
b4.0=           -4.9957 
b4.1=            2.0045 
b5.0=           -1.3523 
b5.1=            0.3713 
b5.2=            0.1693 
b6.0=            1.8056 
b6.1=            0.2224 
b6.2=            0.2265 
b7.0=            3.9801 
b7.1=           -0.0031 
b8.3=            0.5533 
b8.6=           -0.0058 
b9.1=            0.0404 
b9.2=            1.237 
psi.typha=       1.912 

### Prediction equations NOTE THESE SHOULD HAVE LINEAR RELATIONSHIPS IN THESE UNITS
buffp.hat <- b2.0 + b2.1*surr                     # in logits
buffp.hat.pr <- (1/(1+1/(exp(buffp.hat))))*3      # transform from logits to proportions and rescale to (0:3)
hydp.hat <- b3.0 + b3.1*surr
hydp.hat.pr <- (1/(1+1/(exp(hydp.hat))))*6        # transform from logits to proportions and rescale to (0:6)
soilp.hat <- b4.0 + b4.1*buff                     # in logits
soilp.hat.pr <- (1/(1+1/(exp(soilp.hat))))        # transform from logits to proportions and rescale to (0:1)
wetp.hat <- b5.0 + b5.1*hyd + b5.2*surr           # in logits
wetp.hat.pr <- (1/(1+1/(exp(wetp.hat))))*365      # transform from logits to proportions and rescale to (0:1)
cond.hat <- b6.0 + b6.1*soil + b6.2*surr          # in raw units
natr.hat <- b7.0 + b7.1*wet                       # in log units
natr.hat.tr <- exp(natr.hat)                      # transform back from log - poisson model
sphagp.hat <- b8.3*cond +b8.6*wet               # in raw units 
sphagp.hat.pr <- (1/(1+1/(exp(sphagp.hat))))*100  # transform from logits to proportions and rescale to (0:100)
typha.hat <- b9.1*cond + b9.2*((cond-psi.typha)>0)*(cond-psi.typha) 

### Residuals (or raw values if variable is exogenous). (Note, for some reason hard to run these lines in Tinn-R; best to copy & paste.)
surr.res <- surr		  #var 1  note: since surr is not affected by other variables in the model, its values are its residuals
buff.res <- buff-buffp.hat.pr 	  #var 2
hyd.res <- hyd-hydp.hat		  #var 3
soil.res <- soil-soilp.hat.pr 	  #var 4
wet.res <- wet-wetp.hat.pr 	  #var 5
cond.res <- cond-cond.hat 	  #var 6
natr.res <- natr-natr.hat.tr 	  #var 7
sphag.res <- sphag-sphagp.hat.pr  #var 8
typha.res <- typha-typha.hat	  #var 9

# Output residuals to a file if you like.
resid.mat <- data.frame(surr.res,buff.res,hyd.res,soil.res,cond.res,wet.res,natr.res,sphag.res,typha.res)
write.table(resid.mat, file="G:/Documents/data/AcadiaNP/SEMwork/SEM-data-May2011/resid.mat.txt", 
 sep=",", na="NA", quote=FALSE, append=FALSE, row.names=FALSE)

### Computing R-squares for the model.
summary(lm(buff ~buffp.hat.pr))
summary(lm(hyd~hydp.hat.pr))
summary(lm(soil~soilp.hat.pr))
summary(lm(cond~cond.hat))
summary(lm(wet~wetp.hat.pr))
summary(lm(natr~natr.hat.tr))
summary(lm(sphag~sphagp.hat.pr))
summary(lm(typha~typha.hat))

###################################################################
# EXAMINING RESIDUAL RELATIONSHIPS FOR EVIDENCE OF MISSING PATHWAYS
###################################################################
### Step 1: execute function "miss.lnks" so it is in memory. 
# Function "miss.lnks" (courtesy of Donald Schoolmaster Jr.).
miss.lnks<-function(mod){
from.to<-matrix(NA,length(mod),2)
colnames(from.to)<-c("from","to")
for(i in 1:nrow(from.to))from.to[i,]<-strsplit(mod[i],"->")[[1]]
var.list<-NULL
var.list[1]<-from.to[1,1]
k=2
for(i in 2:length(from.to))
if(!from.to[i]%in%var.list){
var.list[k]<-from.to[i]
k=k+1}
all.comb<-outer(var.list,var.list,paste,sep="->")
miss.links<-all.comb[upper.tri(all.comb)][!all.comb[upper.tri(all.comb)]%in%mod&!all.comb[lower.tri(all.comb)]%in%mod]
miss.links<-gsub(">","",miss.links)
return(miss.links)}

### Step 2: Declare a model (this is the initial SE model with ns paths removed).
mod <- c("surr->buff", "surr->hyd", "surr->wet", "surr->cond",
         "buff->soil",
         "hyd->wet", 
         "soil->cond",
         "cond->typha", "cond->sphag",
         "wet->sphag", "wet->natr")
       
### Step 3: Find missing linkages.
missing.1<-miss.lnks(mod)
print("number of missing linkages = "); print(length(missing.1))
print(missing.1)

# List of missing linkages (copied from R console):
#[1] "number of missing linkages = "
#[1] 25
#> print(missing.1)
# [1] "buff-hyd"    "surr-soil"   "hyd-soil"    "buff-cond"   "hyd-cond"    "buff-wet"    "soil-wet"    "cond-wet"    "surr-typha" 
#[10] "buff-typha"  "hyd-typha"   "soil-typha"  "wet-typha"   "surr-sphag"  "buff-sphag"  "hyd-sphag"   "soil-sphag"  "typha-sphag"
#[19] "surr-natr"   "buff-natr"   "hyd-natr"    "soil-natr"   "cond-natr"   "typha-natr"  "sphag-natr" 
>
# Quick examination of Pearson correlations among all residuals
print(cor(resid.mat),digits=2)

# Quick examination of Spearman correlations among all residuals
print(cor(resid.mat,method="spearman"),digits=2)

# Step 4: Examine residual distributions for all endogenous variables
par(mfrow=c(2,2)); hist(buff.res); hist(hyd.res); hist(soil.res); hist(cond.res)
par(mfrow=c(2,2)); hist(wet.res); hist(natr.res); hist(sphag.res); hist(typha.res)

### Step 5: Examine all residual relationships
### Note - for  relationships examined, perform more complete tests as needed (ultimate test is significant coefficient in the model).
### "buff-hyd"
plot(hyd.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(hyd.res, buff.res)  
### "surr-soil"   
plot(soil.res ~ surr.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(soil.res, surr.res)
### "hyd-soil" 
plot(soil.res ~ hyd.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(soil.res, hyd.res)
### buff-cond"   
plot(cond.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(cond.res, buff.res)
### "hyd-cond" --- quadratic? -> no  
plot(cond.res ~ hyd.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(cond.res, hyd.res)
### "buff-wet"    
plot(wet.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(wet.res, buff.res)
### "soil-wet"    
plot(wet.res ~ soil.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(wet.res, soil.res)
### "cond-wet"    
plot(cond.res ~ wet.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(cond.res, wet.res)
### "surr-typha" 
plot(typha.res ~ surr.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res, surr.res)
### "buff-typha" 
plot(typha.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res, buff.res)
### "hyd-typha"  
plot(typha.res ~ hyd.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res, hyd.res)
### "soil-typha" 
plot(typha.res ~ soil.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res, soil.res)
### "wet-typha"  
plot(typha.res ~ wet.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res, wet.res)
### "surr-sphag" --- strong negative relationship
plot(sphag.res ~ surr.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res, surr.res)
### "buff-sphag" 
plot(sphag.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res, buff.res)
### "hyd-sphag"  
plot(sphag.res ~ hyd.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res, hyd.res)
### "soil-sphag" 
plot(sphag.res ~ soil.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res, soil.res)
### "typha-sphag"
plot(sphag.res ~ typha.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res, typha.res)
### "surr-natr"  --- strong negative relationship 
plot(natr.res ~ surr.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, surr.res)
### "buff-natr"  
plot(natr.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, buff.res)
### "hyd-natr"    
plot(natr.res ~ hyd.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, hyd.res)
### "soil-natr"   
plot(natr.res ~ soil.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, soil.res)
### "cond-natr"   
plot(natr.res ~ cond.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, cond.res)
### "typha-natr" --- p=0.07 
plot(natr.res ~ typha.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, typha.res)
### "sphag-natr"  --- clear pos relationship
plot(natr.res ~ sphag.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, sphag.res)


#################################################################
# A DEMONSTRATION PLOT FOR INCLUDING IN THE ARTICLE
#################################################################
# Perform simple tests for residual association. - Use other tests if indicated by plots.
cor.test(soil.res,surr.res) 
cor.test(sphag.res,surr.res)
cor.test(typha.res,surr.res)
cor.test(natr.res,surr.res)

# Replot graphs with lines for those with significant assocations. 
par(mfrow=c(2,2)); 
plot(jitter(soil.res)~jitter(surr.res),pch=19,cex.lab=1.4,cex.axis=1.3,xlab="Land Use",ylab="Soil Disturbance") 
plot(sphag.res~jitter(surr.res),pch=19,cex.lab=1.4,cex.axis=1.3,xlab="Land Use",ylab="Sphagnum"); abline(lm(sphag.res~surr.res)) 
plot(typha.res~jitter(surr.res),pch=19,cex.lab=1.4,cex.axis=1.3,xlab="Land Use",ylab="Typha") 
plot(natr.res~jitter(surr.res),pch=19,cex.lab=1.4,cex.axis=1.3,xlab="Land Use",ylab="Native Richness"); abline(lm(natr.res~surr.res)) 

####################################################################


#############################################################################
# FINAL SE MODEL
#############################################################################
# Creating winbugs program and saving as external file.
sink("pathsimple.final.txt")		
cat("
MODEL LR1 {
	# Likelihoods
	for(i in 1:N) {
	  buff[i] ~ dbin(buffp.hat[i],n.buff)				
		   logit(buffp.hat[i]) <- b2.0 + b2.1*surr[i]
    	  hyd[i] ~ dbin(hydp.hat[i],n.hyd)                      
		   logit(hydp.hat[i]) <- b3.0 + b3.1*surr[i]
	  soil[i] ~ dbern(soilp.hat[i])					
		   logit(soilp.hat[i]) <- b4.0 + b4.1*buff[i]
	  wet[i] ~ dbin(wetp.hat[i],n.wet)				
		   logit(wetp.hat[i]) <- b5.0 +b5.1*hyd[i] +b5.2*surr[i]
	  cond[i] ~ dnorm(cond.hat[i],tau.cond)			
		   cond.hat[i] <- b6.0 + b6.1*soil[i] + b6.2*surr[i]
	  natr[i] ~ dpois(natr.hat[i])				
		   log(natr.hat[i]) <- b7.0 +b7.1*wet[i] +b7.2*surr[i] 
	  sphag[i] ~ dbin(sphagp.hat[i],n.sphag)			
		   logit(sphagp.hat[i]) <- b8.0 +b8.3*cond[i] +b8.6*wet[i] +b8.7*hyd[i] +b8.9*soil[i]
	  typha[i] ~ dnorm(typha.hat[i],tau.typha)              
		   typha.hat[i] <- b9.1*cond[i] + b9.2*step(cond[i]-psi.typha)*(cond[i]-psi.typha)
	}
	# Priors
	b2.0 ~ dnorm(0,0.00001); b2.1 ~ dnorm(0,0.00001) 
	b3.0 ~ dnorm(0,0.00001); b3.1 ~ dnorm(0,0.00001); 
	b4.0 ~ dnorm(0,0.00001); b4.1 ~ dnorm(0,0.00001);
	b5.0 ~ dnorm(0,0.00001); b5.1 ~ dnorm(0,0.00001); b5.2 ~ dnorm(0,0.00001);
	b6.0 ~ dnorm(0,0.00001); b6.1 ~ dnorm(0,0.00001); b6.2 ~ dnorm(0,0.00001); 
  	tau.cond ~ dgamma(0.001,0.001)
	b7.0 ~ dnorm(0,0.00001); b7.1 ~ dnorm(0,0.00001); b7.2 ~ dnorm(0,0.00001);
	b8.0 ~ dnorm(0,0.00001); b8.3 ~ dnorm(0,0.00001); b8.6 ~ dnorm(0,0.00001); 
	b8.7 ~ dnorm(0,0.00001); b8.9 ~ dnorm(0,0.00001);
	b9.1 ~ dnorm(0,0.00001); b9.2 ~ dnorm(0,0.00001); 
	tau.typha ~ dgamma(0.001,0.001); psi.typha ~ dunif(1.5,2.5);

	# Derived Quantities
	prob.9.2.GT0 <- step(b9.2);
	}
",fill=TRUE)
sink()
# end of winbugs code creation

# Creating objects to hand to winbugs.
N=37; n.buff=3; n.wet=365; n.sphag=100; n.hyd=6; 
data = list("N","n.buff","n.wet","n.sphag","n.hyd",
           "surr","buff","soil","hyd","wet","cond","natr","sphag","typha")   
parameters <- c("b2.0","b2.1","b3.0","b3.1","b4.0","b4.1",                          
                "b5.0","b5.1","b5.2",
                "b6.0","b6.1","b6.2",
                "b7.0","b7.1","b7.2",
                "b8.0","b8.3","b8.6","b8.7","b8.9",
                "b9.1","b9.2","psi.typha",
                "prob.9.2.GT0")

inits <- function(){list(b2.0=0,b2.1=0,b3.0=0,b3.1=0,b4.0=0,b4.1=0,                                
                 b5.0=0,b5.1=0,b5.2=0,
                 b6.0=0,b6.1=0,b6.2=0,tau.cond=100,
                 b7.0=0,b7.1=0,b7.2=0,
                 b8.0=0,b8.3=0,b8.6=0,b8.7=0,b8.9=0,
                 b9.1=0,b9.2=1,
                 tau.typha=100,psi.typha=2)}

# MCMC settings to hand to winbugs
 nc <- 3      # number of chains
 ni <- 30000  # number of samples for each chain
 nb <- 1000   # number of samples to discard for burn in
 nt <- 1      # thinning rate

# Hand off to winbugs using "bugs" command. Note: because debug=TRUE, you will need to close winbugs before returning to R.
pathsimple.final.out <- bugs(data=data, inits=inits, parameters.to.save=parameters, model.file="pathsimple.final.txt", 
 n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE, working.directory=getwd())

# Print some basic results
print(pathsimple.final.out,digits=4)  # print parameters

# Plot posterior distributions of interest.
hist(pathsimple.final.out$sims.list$b6.1, col="blue")                  # illustration of histogram of parameter ests
densityplot(pathsimple.final.out$sims.list$b6.1,type="density",lwd=2)  # illustration of density plot


###########################################################
# PREDICTION EQNS FOR FINAL MODEL
###########################################################

### parameter estimates
b2.0=            -3.1022 
b2.1=             2.4396 
b3.0=            -1.5981 
b3.1=             1.1392 
b4.0=            -5.0034 
b4.1=             2.0065 
b5.0=            -1.3522 
b5.1=             0.3711 
b5.2=             0.1696 
b6.0=             1.8053 
b6.1=             0.2220 
b6.2=             0.2267 
b7.0=             3.9931 
b7.1=            -0.0024 
b7.2=            -0.1641 
b8.0=             3.6455 
b8.3=            -1.1331 
b8.6=            -0.0044 
b8.7=            -0.0831 
b8.9=            -0.5694 
b9.1=             0.0452 
b9.2=             1.3543 
psi.typha=        1.9014 

### Prediction equations NOTE THESE SHOULD HAVE LINEAR RELATIONSHIPS IN THESE UNITS
buffp.hat <- b2.0 + b2.1*surr                 # in logits
buffp.hat.pr <- (1/(1+1/(exp(buffp.hat))))*3  # transform from logits to proportions and rescale to (0:3)
hydp.hat <- b3.0 + b3.1*surr
hydp.hat.pr <- (1/(1+1/(exp(hydp.hat))))*6    # transform from logits to proportions and rescale to (0:6)
soilp.hat <- b4.0 + b4.1*buff                 # in logits
soilp.hat.pr <- (1/(1+1/(exp(soilp.hat))))    # transform from logits to proportions and rescale to (0:1)
wetp.hat <- b5.0 + b5.1*hyd + b5.2*surr       # in logits
wetp.hat.pr <- (1/(1+1/(exp(wetp.hat))))*365  # transform from logits to proportions and rescale to (0:1)
cond.hat <- b6.0 + b6.1*soil + b6.2*surr      # in raw log10 units
natr.hat <- b7.0 + b7.1*wet +b7.2*surr        # in natural log units
natr.hat.tr <- exp(natr.hat)                  # transform back from log - poisson model
sphagp.hat <- b8.0 +b8.3*cond +b8.6*wet +b8.7*hyd +b8.9*soil   # in raw units 
sphagp.hat.pr <- (1/(1+1/(exp(sphagp.hat))))*100    # transform from logits to proportions and rescale to (0:100)
typha.hat <- b9.1*cond + b9.2*((cond-psi.typha)>0)*(cond-psi.typha)

### Residuals (or equivalent if exogenous)
surr.res <- surr		  #var 1
buff.res <- buff-buffp.hat.pr 	  #var 2
hyd.res <- hyd-hydp.hat.pr	  #var 3
soil.res <- soil-soilp.hat.pr 	  #var 4
wet.res <- wet-wetp.hat.pr 	  #var 5
cond.res <- cond-cond.hat 	  #var 6
natr.res <- natr-natr.hat.tr 	  #var 7
sphag.res <- sphag-sphagp.hat.pr  #var 8
typha.res <- typha-typha.hat	  #var 9

### Capture residuals in a file
final.res <- data.frame(surr.res,buff.res,hyd.res,soil.res,wet.res,cond.res,natr.res,sphag.res,typha.res)

### R-squares
summary(lm(buff ~buffp.hat.pr))
summary(lm(hyd~hydp.hat.pr))
summary(lm(soil~soilp.hat.pr))
summary(lm(wet~wetp.hat.pr))
summary(lm(cond~cond.hat))
summary(lm(natr~natr.hat.tr))
summary(lm(sphag~sphagp.hat.pr))
summary(lm(typha~typha.hat))

### Examination of predicted versus observed
par(mfrow=c(2,2));plot(buff ~buffp.hat.pr);plot(hyd~hydp.hat.pr);plot(soil~soilp.hat.pr);plot(cond~cond.hat)
par(mfrow=c(2,2));plot(wet~wetp.hat.pr);plot(natr~natr.hat.tr);plot(sphag~sphagp.hat.pr);plot(typha~typha.hat)

### Plotting Typha ~ Cond. for Publication
plot(typha~cond,pch=19,xlab="Log10 Surface Water Conductivity",ylab="Log10 Typha Cover",cex.lab=1.4,cex.axis=1.3)
b9.1= 0.0452; b9.2= 1.3543; psi.typha= 1.9014 #coefs
# if y0=0 when x0=0, y1=? when x1=1.4 and y2=? when x2=1.9 and y3=? when x3=2.7
y1 = b9.1*1.4; print(y1) # y1=0.06325
y2 = b9.1*1.9; print(y2) # y1=0.0859
y3 = y2 + b9.2*(2.7-psi.typha); print(y3) #1.1674
lines(c(1.4,1.9),c(0.063,0.086),type="l",lwd=2) #first segment
lines(c(1.9,2.7),c(0.086,1.1674),type="l",lwd=2) #second segment


######################################################################
# EXAMINING RESIDUAL RELATIONSHIPS FOR MISSING PATHWAYS IN FINAL MODEL
######################################################################
### Step 1: execute function "miss.lnks" so it is in memory. 
# Function "miss.lnks" (courtesy of Donald Schoolmaster Jr.).
miss.lnks<-function(mod){
from.to<-matrix(NA,length(mod),2)
colnames(from.to)<-c("from","to")
for(i in 1:nrow(from.to))from.to[i,]<-strsplit(mod[i],"->")[[1]]
var.list<-NULL
var.list[1]<-from.to[1,1]
k=2
for(i in 2:length(from.to))
if(!from.to[i]%in%var.list){
var.list[k]<-from.to[i]
k=k+1}
all.comb<-outer(var.list,var.list,paste,sep="->")
miss.links<-all.comb[upper.tri(all.comb)][!all.comb[upper.tri(all.comb)]%in%mod&!all.comb[lower.tri(all.comb)]%in%mod]
miss.links<-gsub(">","",miss.links)
return(miss.links)}

### Step 2: Declare a model (this is the final SE model).
mod.f <- c("surr->buff", "surr->hyd", "surr->wet", "surr->cond", "surr->natr",
           "buff->soil",
           "hyd->wet","hyd->sphag", 
           "soil->cond","hyd->sphag",
           "cond->typha", "cond->sphag",
           "wet->sphag", "wet->natr")
       
### Step 3: Find missing linkages.
missing.f<-miss.lnks(mod.f)
print("number of missing linkages = "); print(length(missing.f))
print(missing.f)

# List of missing linkages (copied from R console):
#> print("number of missing linkages = "); print(length(missing.f))
#[1] "number of missing linkages = "
#[1] 23
#> print(missing.f)
# [1] "buff-hyd"    "surr-soil"   "hyd-soil"    "buff-cond"   "hyd-cond"    "buff-wet"    "soil-wet"    "cond-wet"    "buff-natr"  
#[10] "hyd-natr"    "soil-natr"   "cond-natr"   "surr-sphag"  "buff-sphag"  "soil-sphag"  "natr-sphag"  "surr-typha"  "buff-typha" 
#[19] "hyd-typha"   "soil-typha"  "wet-typha"   "natr-typha"  "sphag-typha"

# Quick examination of Pearson correlations among all residuals
print(cor(final.res),digits=2)

# Quick examination of Spearman correlations among all residuals
print(cor(final.res,method="spearman"),digits=2)

# Step 4: Examine residual distributions for all endogenous variables
par(mfrow=c(2,2)); hist(buff.res); hist(hyd.res); hist(soil.res); hist(cond.res)
par(mfrow=c(2,2)); hist(wet.res); hist(natr.res); hist(sphag.res); hist(typha.res)

### Step 5: Examine all residual relationships
### Note - for  relationships examined, perform more complete tests as needed (ultimate test is significant coefficient in the model).
### "buff-hyd"
plot(hyd.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(hyd.res, buff.res)  
### "surr-soil"   
plot(soil.res ~ surr.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(soil.res, surr.res)
### "hyd-soil" 
plot(soil.res ~ hyd.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(soil.res, hyd.res)
### buff-cond"   
plot(cond.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(cond.res, buff.res)
### "hyd-cond"   
plot(cond.res ~ hyd.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(cond.res, hyd.res)
### "buff-wet"    
plot(wet.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(wet.res, buff.res)
### "soil-wet"    
plot(wet.res ~ soil.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(wet.res, soil.res)
### "cond-wet"    
plot(cond.res ~ wet.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(cond.res, wet.res)
### "buff-natr"  
plot(natr.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, buff.res)
### "hyd-natr"    
plot(natr.res ~ hyd.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, hyd.res)
### "soil-natr"   
plot(natr.res ~ soil.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, soil.res)
### "cond-natr"   
plot(natr.res ~ cond.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, cond.res)
### "surr-sphag"  
plot(sphag.res ~ surr.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res, surr.res)
### "buff-sphag"  
plot(sphag.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res, buff.res)
### "soil-sphag"  
plot(sphag.res ~ soil.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res, soil.res)
### "natr-sphag"  
plot(natr.res ~ sphag.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, sphag.res)
### "surr-typha"  
plot(typha.res ~ surr.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res, surr.res)
### "buff-typha" 
plot(typha.res ~ buff.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res, buff.res)
### "hyd-typha"   
plot(typha.res ~ hyd.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res, hyd.res)
### "soil-typha"  
plot(typha.res ~ soil.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res, soil.res)
### "wet-typha"   
plot(typha.res ~ wet.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res, wet.res)
### "natr-typha"  
plot(natr.res ~ typha.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res, typha.res)
### "sphag-typha"
plot(sphag.res ~ typha.res,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res, typha.res)


#######################################
# SCENARIOS IN ECOSPHERE PAPER
#######################################

### THREE SCENARIOS/PREDICTIONS/QUERIES FOR TYPHA CONTROL BASED ON FINAL SE MODEL PARAMETER ESTIMATES
### SCENARIO #1 = NO INTERVENTION; #2 = PREVENT BUFFER INTRUSION & SOIL DISTURBANCE; #3 = SET CONDUCTIVITY VALUES TO REFERENCE VALUES 
# Scenario #2: baseline conductivity should be the values when soil=0 and surr=0, i.e., the intercept of the equation "b6.0"
# param mean   sd      2.5%      25%       50%       75%       97.5% (sd = sd of the mean; sd of sample = 0.0407*sqrt(37)= 0.2476)
# b6.0  1.8055 0.0407  1.7260    1.7790    1.8050    1.8320    1.8860

# predicted distributions for conductivity
cond.predict.1 <- rnorm(200,mean=2.01,sd=0.35)   # simulated values with moments of whole set
cond.predict.2 <- rnorm(200,mean=1.80,sd=0.25)   # simulated values with intercept and sd form MCMC of cond~ buff + soil
cond.predict.3 <- rnorm(200,mean=1.415,sd=0.16)  # simulated values with moments from undisturbed subset 
par(mfrow=c(1,3))
hist(cond.predict.1,col="red",xlim=c(1,3.5),main="")
hist(cond.predict.2,col="red",xlim=c(1,3.5),main="")
hist(cond.predict.3,col="red",xlim=c(1,3.5),main="",breaks=4)

# predicted distributions for typha FIGURE 12 IN MANUSCRIPT FOR ECOSPHERE
typha.predict.1 <- b9.1*cond.predict.1 + b9.2*((cond.predict.1-1.9)>0)*(cond.predict.1-1.9)
typha.predict.2 <- b9.1*cond.predict.2 + b9.2*((cond.predict.2-1.9)>0)*(cond.predict.2-1.9)
typha.predict.3 <- b9.1*cond.predict.3 + b9.2*((cond.predict.3-1.9)>0)*(cond.predict.3-1.9)
typha.predict.1 <- ifelse(typha.predict.1 < 0,0,typha.predict.1)
typha.predict.2 <- ifelse(typha.predict.2 < 0,0,typha.predict.2)
typha.predict.3 <- ifelse(typha.predict.3 < 0,0,typha.predict.3)

# Plot Predicted Outcomes
par(mfrow=c(3,2))
hist(cond.predict.1,col="blue",xlim=c(1,3.5),main="Scenario 1",cex.lab=1.5,cex.axis=1.4,xlab="Observed Conductivity")
hist(typha.predict.1,col="red",xlim=c(0,2),main="Scenario 1",breaks=6,cex.lab=1.5,cex.axis=1.4,xlab="Observed Typha Cover")
hist(cond.predict.2,col="blue",xlim=c(1,3.5),main="Scenario 2",cex.lab=1.5,cex.axis=1.4,xlab="Predicted Conductivity") 
hist(typha.predict.2,col="red",xlim=c(0,2),main="Scenario 2",breaks=4,cex.lab=1.5,cex.axis=1.4,xlab="Predicted Typha Cover")
hist(cond.predict.3,col="blue",xlim=c(1,3.5),main="Scenario 3",cex.lab=1.5,cex.axis=1.4,xlab="Predicted Conductivity") 
hist(typha.predict.3,col="red",xlim=c(0,2),main="Scenario 3",breaks=1,cex.lab=1.5,cex.axis=1.4,xlab="Predicted Typha Cover")

#########################################
### QUERIES IN LIEU OF STANDARDIZED COEFS
# For a reference, see GRACE AND BOLLEN 2005. BULLETIN OF THE ECOLOGICAL SOCIETY OF AMERICA.
#########################################
# basic variable properties
summary(surr);summary(buff);summary(hyd);summary(soil);summary(wet);summary(cond);summary(natr);summary(sphag);summary(typha)

### FOR BUFF ~ SURR ###############################################################################################
# Background Info
mean(buff)					      # mean of raw data
mean(surr)                # raw mean = 0.730
range(surr)               # range = 0 to 3

# Load Coefs
b2.0= -3.0985 						# coefs from final model
b2.1=  2.4376 						# coefs from final model

# Prediction Equations
buff.predict.0 <- b2.0 + b2.1*surr
buff.predict.1 <- b2.0 + b2.1*mean(surr); 	print(buff.predict.1)	  # Scenario1: surr is at its mean value
buff.predict.2 <- b2.0 + b2.1*0; 			      print(buff.predict.2) 	# Scenario2: surr = 0
buff.predict.3 <- b2.0 + b2.1*3; 			      print(buff.predict.3)	  # Scenario3: surr = 3

# using straight probabilities
buff.predict.1pr.raw <- (1/(1+1/(exp(buff.predict.1))));  		print(buff.predict.1pr.raw)
buff.predict.2pr.raw <- (1/(1+1/(exp(buff.predict.2))));  		print(buff.predict.2pr.raw)
buff.predict.3pr.raw <- (1/(1+1/(exp(buff.predict.3))));  		print(buff.predict.3pr.raw)

### SO, WHEN YOU CHANGE SURR ACROSS ITS RANGE, PROBABILITY OF BUFF=MAX VALUE CHANGES FROM 0.043 TO 0.985, WHICH IS 0.942 OF ITS RANGE.


### FOR HYD ~ SURR ###############################################################################################
b3.0=           -1.5993 
b3.1=            1.1400 

hyd.predict.1 <- b3.0 + b3.1*mean(surr); 		print(hyd.predict.1)	# Scenario1: surr is at its mean value
hyd.predict.2 <- b3.0 + b3.1*0; 			print(hyd.predict.2) 	# Scenario2: surr = 0
hyd.predict.3 <- b3.0 + b3.1*3; 			print(hyd.predict.3)	# Scenario3: surr = 3

# using straight probabilities (rather than in observed units)
hyd.predict.1pr.raw <- (1/(1+1/(exp(hyd.predict.1))));  		print(hyd.predict.1pr.raw)
hyd.predict.2pr.raw <- (1/(1+1/(exp(hyd.predict.2))));  		print(hyd.predict.2pr.raw)
hyd.predict.3pr.raw <- (1/(1+1/(exp(hyd.predict.3))));  		print(hyd.predict.3pr.raw)

### SO, WHEN YOU CHANGE SURR ACROSS ITS RANGE, PROBABILITY OF HYD CHANGES FROM 0.168 TO 0.861 OR 0.693 OF ITS RANGE

### FOR SOIL ~ BUFF ###############################################################################################
b4.0=           -5.0525 
b4.1=            2.0273 

soil.predict.1 <- b4.0 + b4.1*mean(buff); 	print(soil.predict.1)	  # Scenario1: buff is at its mean value
soil.predict.2 <- b4.0 + b4.1*0; 			      print(soil.predict.2) 	# Scenario2: buff = 0
soil.predict.3 <- b4.0 + b4.1*3; 			      print(soil.predict.3)	  # Scenario3: buff = 3

# using straight probabilities
soil.predict.1pr.raw <- 1/(1+1/(exp(soil.predict.1)));  		print(soil.predict.1pr.raw)
soil.predict.2pr.raw <- 1/(1+1/(exp(soil.predict.2)));  		print(soil.predict.2pr.raw)
soil.predict.3pr.raw <- 1/(1+1/(exp(soil.predict.3)));  		print(soil.predict.3pr.raw)

### So, when you change buff across its range, probability of buff changes from 0.006 to 0.737 of its range.

### FOR WET ~ SURR + HYD ###########################################################################################
#logit(wetp.hat[i]) <- b5.0 +b5.1*hyd[i] +b5.2*surr[i]
b5.0=           -1.3525 
b5.1=            0.3714 
b5.2=            0.1693 

wet.predict.1   <- b5.0 + b5.1*0 + b5.2*mean(surr);  	print(wet.predict.1) 	# Scenario1: hyd at min (0) and surr at mean
wet.predict.2   <- b5.0 + b5.1*6 + b5.2*mean(surr);  	print(wet.predict.2)    # Scenario2: hyd at max (6) and surr at mean
wet.predict.3   <- b5.0 + b5.1*mean(hyd) + b5.2*0;  	print(wet.predict.3)   	# Scenario3: hyd at mean and surr at min (0)
wet.predict.4   <- b5.0 + b5.1*mean(hyd) + b5.2*3;  	print(wet.predict.4)    # Scenario4: hyd at mean and surr at max (3)

wet.predict.1pr <- 1/(1+1/(exp(wet.predict.1)));  	print(wet.predict.1pr)
wet.predict.2pr <- 1/(1+1/(exp(wet.predict.2)));  	print(wet.predict.2pr)  	
wet.predict.3pr <- 1/(1+1/(exp(wet.predict.3)));  	print(wet.predict.3pr)  	
wet.predict.4pr <- 1/(1+1/(exp(wet.predict.4)));  	print(wet.predict.4pr)  	

# b5.1.std = 0.731 – 0.226 = 0.505
# b5.2.std = 0.480 – 0.357 = 0.123

### FOR COND ~ SURR + SOIL #########################################################################################
# cond.hat[i] <- b6.0 + b6.1*soil[i] + b6.2*surr[i] 
b6.0 =           1.8055 
b6.1=            0.2225 
b6.2=            0.2266 
cond.range <- max(cond)-min(cond); print(cond.range)

cond.predict.1 <- b6.0 + b6.1*min(soil) + b6.2*mean(surr); print(cond.predict.1)	# Scenario1: soil at min (0) and surr at mean
cond.predict.2 <- b6.0 + b6.1*max(soil) + b6.2*mean(surr); print(cond.predict.2)	# Scenario2: soil at max (1) and surr at mean
cond.predict.3 <- b6.0 + b6.1*mean(soil) + b6.2*min(surr); print(cond.predict.3)	# Scenario3: soil at mean and surr at min (0)
cond.predict.4 <- b6.0 + b6.1*mean(soil) + b6.2*max(surr); print(cond.predict.4)	# Scenario4: soil at mean and surr at max (3)

b6.1.std <- (cond.predict.2-cond.predict.1)*(1/1.32); print(b6.1.std)  #here I wish to * rangeX/rangeY
b6.2.std <- (cond.predict.4-cond.predict.3)*(1/1.32); print(b6.2.std)

### FOR NATR ~ SURR + WET #########################################################################################
# log(natr.hat[i]) <- b7.0 +b7.1*wet[i] +b7.2*surr[i] 
b7.0=            3.9928 
b7.1=           -0.0024 
b7.2=           -0.1646 
range(natr) # 62 - 11= 51

natr.predict.1 <- b7.0 + b7.1*min(wet) +b7.2*mean(surr); print(natr.predict.1)   	# Scenario1: wet at min and surr at mean
natr.predict.2 <- b7.0 + b7.1*max(wet) +b7.2*mean(surr); print(natr.predict.2)   	# Scenario2: wet at max and surr at mean
natr.predict.3 <- b7.0 + b7.1*mean(wet) +b7.2*min(surr); print(natr.predict.3)   	# Scenario3: wet at mean and surr at min
natr.predict.4 <- b7.0 + b7.1*mean(wet) +b7.2*max(surr); print(natr.predict.4)   	# Scenario4: wet at mean and surr at max

b7.1.std <- (exp(natr.predict.2)-exp(natr.predict.1))*(1/51); print(b7.1.std)	
b7.2.std <- (exp(natr.predict.4)-exp(natr.predict.3))*(1/51); print(b7.2.std)

### FOR SPHAG ~ COND + WET + HYD + SOIL ### (NOTE: THIS LABEL WAS CORRECTED FROM ORIGINALLY PUBLISHED VERSION ########
# logit(sphagp.hat[i]) <- b8.0 +b8.3*cond[i] +b8.6*wet[i] +b8.7*hyd[i] +b8.9*soil[i]
b8.0=           3.6484 
b8.3=          -1.1346 
b8.6=          -0.0045 
b8.7=          -0.0826 
b8.9=          -0.5691 
cond.range <- max(cond)-min(cond)
wet.range <- max(wet)-min(wet) 

# Scenario1: cond at min then max, all others at mean
sphag.predict.1 <- b8.0 +b8.3*min(cond) +b8.6*mean(wet) +b8.7*mean(hyd) +b8.9*mean(soil); print(sphag.predict.1)
sphag.predict.2 <- b8.0 +b8.3*max(cond) +b8.6*mean(wet) +b8.7*mean(hyd) +b8.9*mean(soil); print(sphag.predict.2)
sphag.predict.1pr <- 1/(1+1/(exp(sphag.predict.1)));  	print(sphag.predict.1pr)
sphag.predict.2pr <- 1/(1+1/(exp(sphag.predict.2)));  	print(sphag.predict.2pr)  	

# Scenario2: wet at min then max, all others at mean
sphag.predict.3 <- b8.0 +b8.3*mean(cond) +b8.6*min(wet) +b8.7*mean(hyd) +b8.9*mean(soil); print(sphag.predict.3)
sphag.predict.4 <- b8.0 +b8.3*mean(cond) +b8.6*max(wet) +b8.7*mean(hyd) +b8.9*mean(soil); print(sphag.predict.4)
sphag.predict.3pr <- 1/(1+1/(exp(sphag.predict.3)));  	print(sphag.predict.3pr)  	
sphag.predict.4pr <- 1/(1+1/(exp(sphag.predict.4)));  	print(sphag.predict.4pr)  	

# Scenario3: hyd at min then max, all others at mean
sphag.predict.5 <- b8.0 +b8.3*mean(cond) +b8.6*mean(wet) +b8.7*min(hyd) +b8.9*mean(soil); print(sphag.predict.5)
sphag.predict.6 <- b8.0 +b8.3*mean(cond) +b8.6*mean(wet) +b8.7*max(hyd) +b8.9*mean(soil); print(sphag.predict.6)
sphag.predict.5pr <- 1/(1+1/(exp(sphag.predict.5)));  	print(sphag.predict.5pr)  	
sphag.predict.6pr <- 1/(1+1/(exp(sphag.predict.6)));  	print(sphag.predict.6pr)  	

# Scenario4: soil at min then max, all others at mean
sphag.predict.7 <- b8.0 +b8.3*mean(cond) +b8.6*mean(wet) +b8.7*mean(hyd) +b8.9*min(soil); print(sphag.predict.7)
sphag.predict.8 <- b8.0 +b8.3*mean(cond) +b8.6*mean(wet) +b8.7*mean(hyd) +b8.9*max(soil); print(sphag.predict.8)
sphag.predict.7pr <- 1/(1+1/(exp(sphag.predict.7)));  	print(sphag.predict.7pr)  	
sphag.predict.8pr <- 1/(1+1/(exp(sphag.predict.8)));  	print(sphag.predict.8pr)  	

b8.3.std <- (sphag.predict.2pr-sphag.predict.1pr); print(b8.3.std)
b8.6.std <- (sphag.predict.4pr-sphag.predict.3pr); print(b8.6.std)
b8.7.std <- (sphag.predict.6pr-sphag.predict.5pr); print(b8.7.std)
b8.9.std <- (sphag.predict.8pr-sphag.predict.7pr); print(b8.9.std)


### FOR TYPHA ~ COND #########################################################################################
# typha.hat[i] <- b9.0 + b9.1*cond[i] + b9.2*step(cond[i]-psi.typha)*(cond[i]-psi.typha)
b9.0=          -0.1579 
b9.1=           0.0775 
b9.2=           1.2080 
psi.typha=       1.9217 
cond.range=max(cond)-min(cond); print(cond.range)
typha.range=max(typha)-min(typha); print(typha.range)

# Scenario1: cond at min then max
typha.predict.1 <- b9.0 +b9.1*min(cond) + b9.2*((min(cond)-psi.typha)>0)*(min(cond)-psi.typha); print(typha.predict.1)
typha.predict.2 <- b9.0 +b9.1*max(cond) + b9.2*((max(cond)-psi.typha)>0)*(max(cond)-psi.typha); print(typha.predict.2)

b9.1and2.std <- (typha.predict.2-typha.predict.1)*(cond.range/typha.range); print(b9.1and2.std)


########################################################################################
# ILLUSTRATION OF THE SOLUTION OF THIS MODEL USING NON-BAYESIAN METHODS
# NOTE, THIS METHOD ONLY WORKS CURRENTLY FOR MODELS THAT DO NOT INCLUDE LATENT VARIABLES
########################################################################################
##### Buffer intrusion modeled as proportional response to Landscape Use.
buff.mirror <- 3-buff
buff.count <- cbind(buff,buff.mirror)
glm.buff<-glm(buff.count ~ surr,family=binomial)
summary(glm.buff)
# Visualize results from linear model.
buffp <- buff/max(buff)  # express buff as proportion for plotting
surrj<-jitter(surr)	# jitter for plotting purposes only
buffpj<-jitter(buffp)	# jitter for plotting purposes only
plot(buffpj ~ surrj, pch=16)
xv <- seq(0,3,0.1)
lines(xv, predict(glm.buff, list(surr=xv),type="response"),lwd=2)
# Estimate an R2
buffp.hat <- predict(glm.buff, list(surr=surr), type="response")
buffp.hatj <- jitter(buffp.hat)
plot(buffpj ~ buffp.hatj,pch=16) # jittered plot of predicted vs observed.
lm.111 <- lm(buffp ~ buffp.hat)
abline(lm.111)                   # add linear prediction line on predicted vs obs plot
summary(lm.111)

##### Hydrologic alteration modeled as proportional response to Landscape Use.
hyd.mirror <- 6-hyd
hyd.count <- cbind(hyd,hyd.mirror)
glm.hyd<-glm(hyd.count ~ surr,family=binomial)
summary(glm.hyd)
# Visualize results from linear model.
hydp <- hyd/max(hyd)  # express hyd as proportion for plotting
surrj<-jitter(surr)	# jitter for plotting purposes only
hydpj<-jitter(hydp)	# jitter for plotting purposes only
plot(hydpj ~ surrj, pch=16)
xv <- seq(0,3,0.1)
lines(xv, predict(glm.hyd, list(surr=xv),type="response"),lwd=2) # note that a more complex model would fit this slightly better
# Estimate R2
hydp.hat <- predict(glm.hyd, list(surr=surr), type="response")
hydp.hatj <- jitter(hydp.hat)
plot(hydpj ~ hydp.hatj,pch=16) # jittered plot of predicted vs observed.
lm.112 <- lm(hydp ~ hydp.hat)
abline(lm.112)                   # add linear prediction line on predicted vs obs plot
summary(lm.112)

##### Soil Disturbance modeled as bernoulli response to Buffer Intrusion.
glm.soil<-glm(soil ~ buff,family=binomial)
summary(glm.soil)
# Visualize results from linear model.
buffj<-jitter(buff)	# jitter for plotting purposes only
soilpj<-jitter(soil)	# jitter for plotting purposes only
plot(soilpj ~ buffj, pch=16)
xv <- seq(0,3,0.1)
lines(xv, predict(glm.soil, list(buff=xv),type="response"),lwd=2) 
# Estimate R2
soilp.hat <- predict(glm.soil, list(buff=buff), type="response")
soilp.hatj <- jitter(soilp.hat)
plot(soilpj ~ soilp.hatj,pch=16) # jittered plot of predicted vs observed.
lm.113 <- lm(soil ~ soilp.hat)
abline(lm.113)                   # add linear prediction line on predicted vs obs plot
summary(lm.113)

##### Flooding Duration modeled as proportional response to Landscape Use and Hydrologic Alteration.
wet.mirror <- 365-wet
wet.count <- cbind(wet,wet.mirror)
glm.wet<-glm(wet.count ~ surr + hyd,family=binomial)
summary(glm.wet)
# Estimate R2
wetp.hat <- predict(glm.wet, list(surr=surr,hyd=hyd), type="response")
lm.114 <- lm(wet ~ wetp.hat)
plot(wet ~ wetp.hat,pch=16)
abline(lm.114)                   # add linear prediction line on predicted vs obs plot
summary(lm.114)

##### Conductivity modeled as linear Gaussian response to Landscape Use and Soil Disturbance.
lm.cond<-lm(cond ~ surr + soil)
summary(lm.cond)
cond.hat <- predict(lm.cond)

##### Native Richness modeled as poisson response to Landscape Use and Flooding Duration.
glm.natr <- glm(natr ~ surr + wet,family=poisson)
summary(glm.natr)
# Estimate R2
natrp.hat <- predict(glm.natr, list(surr=surr,wet=wet), type="response")
lm.115 <- lm(natr ~ natrp.hat)
plot(natr ~ natrp.hat,pch=16)
abline(lm.115)                   
summary(lm.115)

##### Sphagnum modeled as proportional response to cond + wet + hyd + soil.
sphag.mirror <- 100-sphag
sphag.count <- cbind(sphag,sphag.mirror)
glm.sphag<-glm(sphag.count ~ cond + wet + hyd + soil,family=binomial)
summary(glm.sphag)
# Estimate R2
sphagp.hat <- predict(glm.sphag, list(cond=cond,wet=wet,hyd=hyd,soil=soil), type="response")
lm.116 <- lm(sphag ~ sphagp.hat)
summary(lm.116)
plot(sphag ~ sphagp.hat,pch=16)
abline(lm.116)                   

### Typha changepoint regression on conductivity
##########################################
library("segmented")

plot(typha ~ cond,pch=16)                      # look at plot
lines(lowess(cond,typha),lwd=2,lty="dashed")   # add lowess line

# fit single changepoint at cond = 1.9 -- you can play with this to compare change points
fit.lm <-lm(typha ~ cond)
summary(fit.lm)
fit.seg<-segmented(fit.lm, seg.Z= ~ cond, psi=1.9)
summary(fit.seg)  # U1.age is difference in slopes across breakpoint
slope(fit.seg) # estimates of slopes
# plot predicted vs observed
typha.hat <- predict(fit.seg)
plot(typha ~ typha.hat,pch=16)
lm.117 <- lm(typha ~ typha.hat)
abline(lm.117)
summary(lm.117)

#####################################
# EXAMINING RESIDUALS FROM ML METHODS
#####################################

### Residuals (or equivalent if exogenous)   # for some reason these commands had to be pasted to the console when using Tinn-R
surr.res2  <- surr		    #var 1
buff.res2  <- buff - buffp.hat      #var 2
hyd.res2   <- hydp - hydp.hat	    #var 3
soil.res2  <- soil - soilp.hat      #var 4
wet.res2   <- wet - wetp.hat        #var 5
cond.res2  <- cond - cond.hat	    #var 6
natr.res2  <- natr - natrp.hat 	    #var 7
sphag.res2 <- sphag - sphagp.hat    #var 8
typha.res2 <- typha - typha.hat	    #var 9

### Plotting residuals for conditionally independent pairs (as performed for Final Model in Bayesian section)
### "buff-hyd"
plot(hyd.res2 ~ buff.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(hyd.res2, buff.res2)  
### "surr-soil"   
plot(soil.res2 ~ surr.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(soil.res2, surr.res2)
### "hyd-soil" 
plot(soil.res2 ~ hyd.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(soil.res2, hyd.res2)
### buff-cond"   
plot(cond.res2 ~ buff.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(cond.res2, buff.res2)
### "hyd-cond"   
plot(cond.res2 ~ hyd.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(cond.res2, hyd.res2)
### "buff-wet"    
plot(wet.res2 ~ buff.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(wet.res2, buff.res2)
### "soil-wet"    
plot(wet.res2 ~ soil.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(wet.res2, soil.res2)
### "cond-wet"    
plot(cond.res2 ~ wet.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(cond.res2, wet.res2)
### "buff-natr"  
plot(natr.res2 ~ buff.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res2, buff.res2)
### "hyd-natr"    
plot(natr.res2 ~ hyd.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res2, hyd.res2)
### "soil-natr"   
plot(natr.res2 ~ soil.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res2, soil.res2)
### "cond-natr"   
plot(natr.res2 ~ cond.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res2, cond.res2)
### "surr-sphag" # indicates possible missing linkage; however, surr does not stay in the model for sphag, so link not justified
plot(sphag.res2 ~ surr.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res2, surr.res2)
### "buff-sphag"   
plot(sphag.res2 ~ buff.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res2, buff.res2)
### "soil-sphag"  
plot(sphag.res2 ~ soil.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res2, soil.res2)
### "natr-sphag"  
plot(natr.res2 ~ sphag.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res2, sphag.res2)
### "surr-typha"  
plot(typha.res2 ~ surr.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res2, surr.res2)
### "buff-typha" 
plot(typha.res2 ~ buff.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res2, buff.res2)
### "hyd-typha"   
plot(typha.res2 ~ hyd.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res2, hyd.res2)
### "soil-typha"  
plot(typha.res2 ~ soil.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res2, soil.res2)
### "wet-typha"   
plot(typha.res2 ~ wet.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(typha.res2, wet.res2)
### "natr-typha"  
plot(natr.res2 ~ typha.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(natr.res2, typha.res2)
### "sphag-typha"
plot(sphag.res2 ~ typha.res2,pch=19,cex.lab=1.4,cex.axis=1.3)
cor.test(sphag.res2, typha.res2)
