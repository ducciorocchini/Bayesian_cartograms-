# Load packages
library(rgdal)
library(lme4)
library(maptools)
library(R2jags)
library(pROC)

# Set working directory
setwd('./')

# Load county spatial data
county<-readOGR('USA_counties_4326.sqlite')

# Load L.japonica presence per county
pa<-read.csv('bioclim_table.csv',header=TRUE)

# Merge spatial data with observation and bioclim variables
dat.pa<-merge(county,pa,by=c("name","state_name"))

# Remove counties without bioclim data
dat.pa<-subset(dat.pa,!is.na(mean_bio6))

# Subset the dataset to 10% of representated states
set.seed(26)
sub.pa<-subset(dat.pa,state_name%in%sample(dat.pa$state_name, 10),replacement=FALSE)

# Rescale bio6 and center bio12
sub.pa@data[,c("mean_bio6")]<-sub.pa@data[,c("mean_bio6")]/10
sub.pa@data[,c("mean_bio12")]<-sub.pa@data[,c("mean_bio12")]-1000
dat.pa@data[,c("mean_bio6")]<-dat.pa@data[,c("mean_bio6")]/10
dat.pa@data[,c("mean_bio12")]<-dat.pa@data[,c("mean_bio12")]-1000

# Add bio6^2
sub.pa@data$mean_bio6.2<-(sub.pa@data[,c("mean_bio6")]-5)^2
dat.pa@data$mean_bio6.2<-(dat.pa@data[,c("mean_bio6")]-5)^2

####LMM#####
# Run lmm with state as group factor
lmer_jap<-glmer(japhoney~scale(mean_bio6.2,center=FALSE)+scale(mean_bio12,center=FALSE)+(1|state_name),family="binomial",data=sub.pa@data)

# Prepare dataset for prediction
dat.pred<-cbind.data.frame(state_name=dat.pa@data$state_name,mean_bio6.2=scale(dat.pa@data$mean_bio6.2,center=FALSE),mean_bio12=scale(dat.pa@data$mean_bio12,center=FALSE))

# Predict with lmm in the scale of the linear predictor
lmer_jap_pred<-predict(lmer_jap,newdata=dat.pred,type = "response",allow.new.levels = TRUE)

# Calculate 95% CI of the linear predictor
predFun <- function(fit) {
  predict(fit,dat.pred,type = "response",allow.new.levels = TRUE)
}

bb <- bootMer(lmer_jap,nsim=99,parallel="multicore",FUN=predFun, ncpus=8,use.u = FALSE,.progress="txt")
CI95<-apply(bb$t,2,quantile,p=c(0.025,0.975),na.rm=TRUE)
se95<-(CI95[2,]-CI95[1,])

# Add LMM prediction and standard error to the map
dat.pa@data$lmer_p <- lmer_jap_pred
dat.pa@data$lmer_se <- se95

####BHLMM#####
# Define the model in JAGS language
cat("model {

# Priors
mu_int~dnorm(0, tau_int) # Mean hyperparameter for random intercepts
sigma_int~dlnorm(0,1) # SD hyperparameter for random intercepts
tau_int <- pow(sigma_int,-2)

beta_bio62~dnorm(-0.5, 1/pow(0.4,2)) # Common slope
beta_bio12~dnorm(0.1, 1/pow(0.1,2)) # Common slope

# Likelihood
for (i in 1:n_obs) {
    y[i]~dbern(mu[i]) # The actual (random) responses
    logit(mu[i]) <- alpha[states[i]]+ beta_bio62*bio62[i] + beta_bio12*bio12[i] # Expectation
}

for(i in 1:n_states) {
alpha[i]~dnorm(mu_int, tau_int) # Random intercepts

}

}", fill=TRUE, file="bhlmm_japhoney.txt")

# Bundle data for BHLMM testing and training
jags_data <- list(y=as.numeric(dat.pa$japhoney), states=as.factor(dat.pa$state_name), bio62=as.numeric(scale(dat.pa$mean_bio6.2,center=FALSE)), bio12=as.numeric(scale(dat.pa$mean_bio12,center=FALSE)), n_states=length(unique(dat.pa$state_name)), n_obs=length(dat.pa$japhoney))

#Set NA for values not in sub.pa
jags_data[[1]][which(1:length(jags_data[[1]])%not in%as.integer(row.names(sub.pa@data)))]<-NA

# Drop unused levels (JAGS requires continous numeration of levels)
jags_data$states <- droplevels(jags_data$states)
jags_data$states <- as.integer(jags_data$states)

# Set parameters to estimate
params <- c("alpha","mu_int","sigma_int","beta_bio62","beta_bio12")
params_pred <- c("mu")

# MCMC settings
ni <- 10000; nb <- 2000; nt <- 5; nc <- 3

# Estimate model parameters
outj <- jags.parallel(data=jags_data,model.file="bhlmm_japhoney.txt", n.chains=nc,parameters.to.save = params,n.iter=ni,export_obj_names=c("ni","nc"),digits=3)

# Predict
outjp <- jags.parallel(data=jags_data,model.file="bhlmm_japhoney.txt", n.chains=nc,parameters.to.save = params_pred,n.iter=ni,export_obj_names=c("ni","nc"),digits=3,DIC=FALSE)

# Add BHLMM prediction to the spatial data frame
dat.pa@data$blmm_p <-  as.numeric(outjp[[2]]$mean[[1]])
dat.pa@data$blmm_sd <- as.numeric(outjp[[2]]$sd[[1]])

# Build ROCs for LMM and BHLMM and plot them 
g <- roc(japhoney ~ lmer_p, data = dat.pa)
g1 <- roc(japhoney ~ blmm_p, data = dat.pa)

g2 <- ggroc(list(LMM=g, BHLMM=g1), linetype=1, legacy.axes=TRUE) + scale_colour_discrete("Model") + theme(legend.direction = "horizontal", legend.position = "top",legend.title=element_blank()) + geom_abline(slope=1,intercept=0,linetype=2)