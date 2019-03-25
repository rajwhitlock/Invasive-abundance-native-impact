#Meta-analysis code and analysis framework written by Raj Whitlock (RW) 08/01/2017
#Code last updated by Bethany Bradley and RW 03/24/2019
#Overview: this code allows raw data pre-processing to calculate "slopes" effect sizes, meta-analysis of processed effect-size data and figure creation from model output
#You will need raw input files "AvI_data_sheet.txt" and "AvI_attribute_sheet.txt" as well as the R script 'Bar_plot_error_bars.R'. Place the datasets in your working directory
##################################################################
##################################################################
##################################################################
##################################################################

source("/Directory/path/to/script/Bar_plot_error_bars.R")

setwd("/Directory/path/to/working_directory")

#load necessary libraries
library(MASS)
library(MCMCglmm)
library(metafor)

# SECTION 1. The data
full_data <- read.delim("AvI_data_sheet.txt")
aie_db <- read.delim("AvI_attribute_sheet.txt", header=T)



##################################################################
##################################################################
##################################################################
# SECTION 2. Needed utility functions:
# This is a function to calculate approximate Bayesian p-vales from posterior samples from MCMCglmm models, following Hadfield JD (2010) MCMC methods for multi-response generalized linear mixed models: the MCMCglmm R package. Journal of Statistical Software 33, 1–22.

pMCMC <- function(data){
res <- 2*(min(length(which(data >0)),length(which(data <0))))/length(data)
if (res==0){res <- "<0.001"}
return(as.character(res)) 
}



##################################################################
##################################################################
##################################################################

# SECTION 3. Data pre-processing to calculate effect sizes


# First, define some "columns to pass" (c2p) through data pre-processing that are cut out and joined back into the effect size dataset prior to analysis
c2p <- c(which(names(data)=="Article_ID"),which(names(data)=="Study_ID"),which(names(data)=="Art_stud"))

##################################################################
##################################################################
####################START OF FUNCTION CREATING "SLOPES" EFFECT SIZES########

# Function to process the raw data and extract regression coefficients and their standard error (effect sizes). The argument columns.to.pass specifies article-level variates that should be retained, excluding the unique article identifier, which will be retained automatically

avi.es.calc2 <- function (data,columns.to.pass, rsc.x = T, rsc.y = T){
# data = input data
# columns.to.pass are study-level variates embedded in the data that we want to keep (they are cut out, and added back in after effect size calculation)
# rsc.x, and rsc.y are options to rescale x and y to c(0, 1), default = T
# the raw x data are always centred

avi.cov <- data[,c2p]
avi.cov2 <- split(avi.cov,f = list(avi.cov[,1],avi.cov[,2]),drop=T)
avi.cov3 <- lapply(avi.cov2,function(x){
	x[1,]
})
avi.cov <- do.call(rbind,avi.cov3)

xx <- split(data,list(data$Article_ID,data$Study_ID),drop=T)	
yy <- lapply(xx, function(y){
#rescaling the data between 0-1 (if rsc.x and rsc.y are True)
	n <- dim(y)[1]
	if (rsc.x == T){
		if(sign (min (y$Abundance_Invader))==-1){offset.x <- min (y$Abundance_Invader)} else {offset.x <- 0}
		inv.recentred <- (y$Abundance_Invader - offset.x)/(max(y$Abundance_Invader) - offset.x)
		inv.recentred <- inv.recentred - mean(inv.recentred)
	} else {
		inv.recentred <- y$Abundance_Invader - mean(y$Abundance_Invader)
	}
	
	if (rsc.y == T){
		if(sign (min (y$Response))==-1){offset.y <- min (y$Response)} else {offset.y <- 0}
		y.resp <- (y$Response - offset.y) / (max(y$Response) - offset.y)
	} else {y.resp <- y$Response}
	
	#mev <- 1/(n-3)
	
	# if a study only has 3 effective invasive data points, we can't compute the polynomial slope, so I give the option to estimate linear only for these ones
	if (n > 3 && length(unique(y$Abundance_Invader)) > 3){
		
	m1 <- lm(y.resp ~ poly(inv.recentred,2, raw=T)) #second order polynomial fit between response and invasive
	# Raw polynomials were chosen so that a consistently scaled model matrix applied in each of the studies, allowing comparability and synthesis of effect sizes from different studies
	
	a.1 <- summary(m1)$coef[1,1] #estimate
	a.2 <- summary(m1)$coef[1,2] #St Error
	
	a.3 <- summary(m1)$coef[2,1] #estimate
	a.4 <- summary(m1)$coef[2,2] #St Error
	
	a.5 <- summary(m1)$coef[3,1] #estimate
	a.6 <- summary(m1)$coef[3,2] #St Error
	c.m1 <- cov2cor(vcov(m1))[3,2]	# calculate within-study correlation of linear and polynomial predictors within studies
	
	}   else    {
	
	m1 <- lm(y.resp ~ inv.recentred)
	a.1 <- summary(m1)$coef[1,1]
	a.2 <- summary(m1)$coef[1,2]
	
	a.3 <- summary(m1)$coef[2,1]
	a.4 <- summary(m1)$coef[2,2]
	
	a.5 <- NA # no polynomial slope estimated
	a.6 <- NA # no polynomial slope estimated
	c.m1 <- NA
	
	}

	return(data.frame("Article_ID" = as.character(y$Article_ID[1]), "Study_ID" = as.character(y$Study_ID[1]), "Int.u" = a.1, "Int.s" = a.2, "Lin.u" = a.3, "Lin.s" = a.4, "Pol.u" = a.5, "Pol.s" = a.6, "cor1.2" = c.m1))
})

yy2 <- do.call(rbind,yy)

res <- merge (avi.cov,yy2)
return (res)
}


##################################################################
####################END OF FUNCTION CREATING SLOPES EFFECT SIZES########
##################################################################

# SECTION 4. Now run the function and create the effect size dataset:

avi.ES.s <- avi.es.calc2(data,columns.to.pass=c2p, rsc.x = T, rsc.y = T)

#database internal join of effect estimates with article level covariates
avi.ES.s <- merge(avi.ES.s, aie_db)

avi.ES.s <- avi.ES.s[-which(is.na(avi.ES.s$Pol.u)),] #removes studies where polynomial slopes could not be estimated (<4 data points)
#Two studies removed: BRAN2010_1 and TERA2011_2

#Two studies with problematic polynomial slopes identified: BRAN2010_3 (only two response values, poly.u = 60)
#TERA2011_4 poly.u = -464 (only five lines in study).  Exclude both of these
avi.ES.s <- avi.ES.s[-which(avi.ES.s$Art_stud=="BRAN2010_3"),]
avi.ES.s <- avi.ES.s[-which(avi.ES.s$Art_stud=="TERA2011_4"),]

##################################################################
##################################################################
##################################################################


# SECTION 5. Meta-analysis models. The code in this section specifies models that can be run on the effect size data, and code for producing figures from model output. Take care, the correspondence between figure numbers indicated in this code and figure numbers in the manuscript and its supplementary materials may need checking

# now for the analysis, random effects model with article level random effects

#We need this prior: uniform prior on the standard deviation of the random effects, for both residual variance (R structure), and study level variance (G structure)
prior <- list(R=list(V = 1e-10, nu = -1), G = list(G1 = list(V = 1e-10, nu = -1)))

##################################################################
##################################################################
## Analysis 1: Global meta-analysis ## (Not included in manuscript)

# random effects meta-analysis with additional random effects for study, no fixed effects bar the intercept

ms1a <- MCMCglmm(Int.u ~ 1, random = ~ Article_ID, mev = avi.ES.s$Int.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES.s)
ms1b <- MCMCglmm(Lin.u ~ 1, random = ~ Article_ID, mev = avi.ES.s$Lin.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES.s)
ms1c <- MCMCglmm(Pol.u ~ 1, random = ~ Article_ID, mev = avi.ES.s$Pol.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES.s)
#mev is a vector of measurement error variance for each data point

# Effect size, direction and statistical significance of linear term
summary(ms1b)

# Effect size, direction and statistical significance of polynomial term
summary (ms1c)

# Plot the functional response curve (and 95% credible zone) from meta-analysis of slopes

#prediction interval on the centred x scale (-0.5, 0.5).  Subsequent plots are on the original rescaled (0, 1 scale) (0.5 is added to predicted values on the x axis)
ix <- (-500:500)/1000

# extract posteriors for the average regression coefficients (intercept, x, x^2)
ms1.sol <- cbind(ms1a$Sol[, "(Intercept)"],ms1b$Sol[, "(Intercept)"],ms1c$Sol[, "(Intercept)"])
# split by row, for lapply convenience later...
ms1.sol <- split(ms1.sol,1:1000)

# matrix to catch results for 95% credible zone
res <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms1.sol, function(y){
    ny <- y[3]*(ix[i]^2) + y[2]*ix[i] + y[1]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res[i,] <- c(pred1.int,pred1.u)
}

######## Proportional change in native responses ###################################
######## over the typical range in invasive species' ###############################
######## abundance investigated in the literature ##################################
# note: changes are on average over meta-analysed studies, the response in individual studies varies.
res[1000,3] - res[1,3]
# On average, there was a 23.2% decrease in native responses over the typical range in invasive species abundance 

##########################################################

# plot the functional response curve and 95% credible zone

## Figure X: Slopes Global #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res[,1],res[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res[,3],type="l", lwd = 2)


##################################################################
##################################################################
## Analysis 2: Community vs. Population level response
##################################################################
#avi.ES.s <- avi.ES.s[-which (avi.ES.s$Study_type == "Spatial"),] #double check to make sure results still hold when spatial studies excluded. They do.


avi.ES2.s <- avi.ES.s[-which(avi.ES.s$Response_type == "Other"),]




ms2a <- MCMCglmm(Int.u ~ Multi_spp_resp, random = ~ Article_ID, mev = avi.ES2.s$Int.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES2.s)
ms2b <- MCMCglmm(Lin.u ~ Multi_spp_resp, random = ~ Article_ID, mev = avi.ES2.s$Lin.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES2.s)
ms2c <- MCMCglmm(Pol.u ~ Multi_spp_resp, random = ~ Article_ID, mev = avi.ES2.s$Pol.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES2.s)

summary(ms2a); summary(ms2b); summary(ms2c)

# Community level intercept term is significantly greater than from zero p < 0.001
# Population level intercept term is significantly lower than community level intercept term p < 0.001

# Community level linear term is significantly less than zero p < 0.001
# Population level linear term does not differ significantly from community level linear term p = 0.282

# Community level polynomial term does not differ significantly from zero p = 0.330
# Population level polynomial term is significantly greater than community level polynomial term p = 0.026


## Plot the functional response curve (and 95% credible zone) ##########################
########################################################################################

# extract posteriors for the average regression coefficients (intercept, x, x^2)
ms2.sol <- cbind(
	ms2a$Sol[, "(Intercept)"] + ms2a$Sol[, "Multi_spp_respSINGLE"],
	ms2b$Sol[, "(Intercept)"] + ms2b$Sol[, "Multi_spp_respSINGLE"],
	ms2c$Sol[, "(Intercept)"] + ms2c$Sol[, "Multi_spp_respSINGLE"],
	ms2a$Sol[, "(Intercept)"],
	ms2b$Sol[, "(Intercept)"],
	ms2c$Sol[, "(Intercept)"])


# Estimate for Population, linear term
mean(mcmc(ms2.sol)[,2])
# Estimate for Population, polynomial term
mean(mcmc(ms2.sol)[,3])
# Estimate for Community, linear term
mean(mcmc(ms2.sol)[,5])
# Estimate for Community, polynomial term
mean(mcmc(ms2.sol)[,6])

# p-value for Population, linear term, comparing to zero
pMCMC(mcmc(ms2.sol)[,2])
# <0.001
# p-value for Population, polynomial term, comparing to zero
pMCMC(mcmc(ms2.sol)[,3])
# =0.002
# p-value for Community, linear term, comparing to zero
pMCMC(mcmc(ms2.sol)[,5])
# <0.001
# p-value for Community, polynomial term, comparing to zero
pMCMC(mcmc(ms2.sol)[,6])
# =0.330


# split by row, for lapply convenience later...
ms2.sol <- split(ms2.sol,1:1000)
ix <- (-500:500)/1000

####### (i) Population level predictions and credible zone ###############################

# matrix to catch results for 95% credible zone
res2a <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms2.sol, function(y){
    ny <- y[3]*(ix[i]^2) + y[2]*ix[i] + y[1]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res2a[i,] <- c(pred1.int,pred1.u)
}

####### (ii) Community level predictions and credible zone ###############################

# matrix to catch results for 95% credible zone
res2b <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms2.sol, function(y){
	# Note indexing within the square brackets refers to correct columns in ms2.sol
    ny <- y[6]*(ix[i]^2) + y[5]*ix[i] + y[4]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res2b[i,] <- c(pred1.int,pred1.u)
}


######## Proportional change in native responses ###################################
######## over the typical range in invasive species' ###############################
######## abundance investigated in the literature ##################################
# note: changes are on average over meta-analysed studies, the response in individual studies varies.

res2a[1000,3] - res2a[1,3]
# There was a 19.9% decrease in population-level native responses over the typical range in invasive species abundance investigated in the literature

res2b[1000,3] - res2b[1,3]
# There was a 24.6% decrease in community-level native responses over the typical range in invasive abundance investigated in the literature

##########################################################

## Figure 2b: Slopes Population #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res2a[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Population response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res2a[,1],res2a[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res2a[,3],type="l", lwd = 2)



## Figure 2d: Slopes Community #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res2b[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Community response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res2b[,1],res2b[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res2b[,3],type="l", lwd = 2)



##################################################################
## Analysis 3: Community vs. Population by trophic level 
##################################################################

avi.ES3.s <- avi.ES2.s[-c(which(avi.ES2.s$Trophic_level=="Above/Intra"),which(avi.ES2.s$Trophic_level=="Below/Intra"),which(avi.ES2.s$Trophic_level=="Mixed"),which(avi.ES2.s$Trophic_level=="Unknown")),]
# remove empty levels (clean up)
avi.ES3.s <- droplevels (avi.ES3.s)


ms4a <- MCMCglmm(Int.u ~ Trophic_level + Multi_spp_resp, random = ~ Article_ID, mev = avi.ES3.s$Int.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3.s)
ms4b <- MCMCglmm(Lin.u ~ Trophic_level + Multi_spp_resp, random = ~ Article_ID, mev = avi.ES3.s$Lin.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3.s)
ms4c <- MCMCglmm(Pol.u ~ Trophic_level + Multi_spp_resp, random = ~ Article_ID, mev = avi.ES3.s$Pol.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3.s)

summary(ms4b)

# The linear term for Trophic category Above is significantly less than zero p < 0.001
# The linear term for Intra is significantly greater than that for Above  p = 0.002
# The linear term for Below is significantly greater than that for Above  p < 0.001

# set reference level of Trophic_level to intra, to examine intra vs. below
avi.ES3.s$Trophic_level <- relevel (avi.ES3.s$Trophic_level, ref = 3)

#model m4b2
m4b2 <- MCMCglmm(Lin.u ~ Trophic_level + Multi_spp_resp, random = ~ Article_ID, mev = avi.ES3.s$Lin.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3.s)
summary(m4b2)

# Linear effect size for intra is significantly lower than (more negative) than that for  below p < 0.001

# VERY IMPORTANT: Reset the Trophic_level variable to original indexing
avi.ES3.s$Trophic_level <- relevel (avi.ES3.s$Trophic_level, ref = 2)

summary(ms4c)
# The polynomial term for Trophic category Above is significantly greater than zero p = 0.002
# The polynomial term for Intra is significantly lower than that for Above  p = 0.018
# The polynomial term for Below is significantly lower than that for Above  p = 0.026


# set reference level of Trophic_level to intra, to examine intra vs. below
avi.ES3.s$Trophic_level <- relevel (avi.ES3.s$Trophic_level, ref = 3)

#model m4c2
m4c2 <- MCMCglmm(Pol.u ~ Trophic_level + Multi_spp_resp, random = ~ Article_ID, mev = avi.ES3.s$Pol.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3.s)
summary(m4c2)

# Polynomial effect size for intra does not differ significantly from that for below p = 0.586

# VERY IMPORTANT: Reset the Trophic_level variable to original indexing
avi.ES3.s$Trophic_level <- relevel (avi.ES3.s$Trophic_level, ref = 2)




## Plot the functional response curve (and 95% credible zone) ##########################
########################################################################################

# extract posteriors for the average regression coefficients (intercept, x, x^2)
ms4.sol <- cbind(
	ms4a$Sol[, "(Intercept)"] + ms4a$Sol[, "Multi_spp_respSINGLE"],
	ms4b$Sol[, "(Intercept)"] + ms4b$Sol[, "Multi_spp_respSINGLE"],
	ms4c$Sol[, "(Intercept)"] + ms4c$Sol[, "Multi_spp_respSINGLE"],
	ms4a$Sol[, "(Intercept)"] + ms4a$Sol[, "Multi_spp_respSINGLE"]+ ms4a$Sol[, "Trophic_levelIntra"],
	ms4b$Sol[, "(Intercept)"] + ms4b$Sol[, "Multi_spp_respSINGLE"]+ ms4b$Sol[, "Trophic_levelIntra"],
	ms4c$Sol[, "(Intercept)"] + ms4c$Sol[, "Multi_spp_respSINGLE"]+ ms4c$Sol[, "Trophic_levelIntra"],
	ms4a$Sol[, "(Intercept)"] + ms4a$Sol[, "Multi_spp_respSINGLE"] + ms4a$Sol[, "Trophic_levelBelow"],
	ms4b$Sol[, "(Intercept)"] + ms4b$Sol[, "Multi_spp_respSINGLE"] + ms4b$Sol[, "Trophic_levelBelow"],
	ms4c$Sol[, "(Intercept)"] + ms4c$Sol[, "Multi_spp_respSINGLE"] + ms4c$Sol[, "Trophic_levelBelow"],
	ms4a$Sol[, "(Intercept)"],
	ms4b$Sol[, "(Intercept)"],
	ms4c$Sol[, "(Intercept)"],
	ms4a$Sol[, "(Intercept)"]+ ms4a$Sol[, "Trophic_levelIntra"],
	ms4b$Sol[, "(Intercept)"]+ ms4b$Sol[, "Trophic_levelIntra"],
	ms4c$Sol[, "(Intercept)"]+ ms4c$Sol[, "Trophic_levelIntra"],
	ms4a$Sol[, "(Intercept)"] + ms4a$Sol[, "Trophic_levelBelow"],
	ms4b$Sol[, "(Intercept)"] + ms4b$Sol[, "Trophic_levelBelow"],
	ms4c$Sol[, "(Intercept)"] + ms4c$Sol[, "Trophic_levelBelow"])

# Estimates for 18 regression term effect sizes
apply(mcmc(ms4.sol), 2, mean)
# Ordering of estimates is as follows, these contain information to add to plots
# 1 Population, Above, Intercept
# 2 Population, Above, Linear
# 3 Population, Above, Polynomial
# 4 Population, Intra, Intercept
# 5 Population, Intra, Linear
# 6 Population, Intra, Polynomial
# 7 Population, Below, Intercept
# 8 Population, Below, Linear
# 9 Population, Below, Polynomial
# 10 Community, Above, Intercept
# 11 Community, Above, Linear
# 12 Community, Above, Polynomial
# 13 Community, Intra, Intercept
# 14 Community, Intra, Linear
# 15 Community, Intra, Polynomial
# 16 Community, Below, Intercept
# 17 Community, Below, Linear
# 18 Community, Below, Polynomial

# p-values for 18 regression term effect sizes, comparing to zero
apply(mcmc(ms4.sol), 2, pMCMC)
# Ordering is as in comments immediately preceding, some of these p-values to be added to plots as star symbols



# split by row, for lapply convenience later...
ms4.sol <- split(ms4.sol,1:1000)
ix <- (-500:500)/1000

####### (i) Population, Above: predictions and credible zone #############################

# matrix to catch results for 95% credible zone
res4a <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms4.sol, function(y){
    ny <- y[3]*(ix[i]^2) + y[2]*ix[i] + y[1]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res4a[i,] <- c(pred1.int,pred1.u)
}

####### (ii) Population, Intra: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res4b <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms4.sol, function(y){
	# Note indexing within the square brackets refers to correct columns in ms4.sol
    ny <- y[6]*(ix[i]^2) + y[5]*ix[i] + y[4]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res4b[i,] <- c(pred1.int,pred1.u)
}

####### (iii) Population, Below: predictions and credible zone ##########################

# matrix to catch results for 95% credible zone
res4c <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms4.sol, function(y){
	# Note indexing within the square brackets refers to correct columns in ms4.sol
    ny <- y[9]*(ix[i]^2) + y[8]*ix[i] + y[7]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res4c[i,] <- c(pred1.int,pred1.u)
}

####### (iv) Community, Above: predictions and credible zone ############################

# matrix to catch results for 95% credible zone
res4d <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms4.sol, function(y){
	# Note indexing within the square brackets refers to correct columns in ms4.sol
    ny <- y[12]*(ix[i]^2) + y[11]*ix[i] + y[10]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res4d[i,] <- c(pred1.int,pred1.u)
}

####### (v) Community, Intra: predictions and credible zone #############################

# matrix to catch results for 95% credible zone
res4e <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms4.sol, function(y){
	# Note indexing within the square brackets refers to correct columns in ms4.sol
    ny <- y[15]*(ix[i]^2) + y[14]*ix[i] + y[13]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res4e[i,] <- c(pred1.int,pred1.u)
}

####### (vi) Community, Below: predictions and credible zone #############################

# matrix to catch results for 95% credible zone
res4f <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms4.sol, function(y){
	# Note indexing within the square brackets refers to correct columns in ms4.sol
    ny <- y[18]*(ix[i]^2) + y[17]*ix[i] + y[16]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res4f[i,] <- c(pred1.int,pred1.u)
}

######## Proportional change in native responses ###################################
######## over the typical range in invasive species' ###############################
######## abundance investigated in the literature ##################################
# note: changes are on average over meta-analysed studies, the response in individual studies varies.

res4a[1000,3] - res4a[1,3]
# Where invasive species occupied a higher trophic level, there was a 44.0% decrease in population-level native responses over the typical range in invasive abundance investigated in the literature
res4b[1000,3] - res4b[1,3]
# Where invasive species occupied the same trophic level, there was a 19.7% decrease in population-level native responses over the typical range in invasive abundance investigated in the literature
res4c[1000,3] - res4c[1,3]
# Where invasive species occupied a lower trophic level, there was a 0.5% increase in population-level native responses over the typical range in invasive abundance investigated in the literature
res4d[1000,3] - res4d[1,3]
# Where invasive species occupied a higher trophic level, there was a 52.0% decrease in community-level native responses over the typical range in invasive abundance investigated in the literature
res4e[1000,3] - res4e[1,3]
# Where invasive species occupied the same trophic level, there was a 27.8% decrease in community-level native responses over the typical range in invasive abundance investigated in the literature
res4f[1000,3] - res4f[1,3]
# Where invasive species occupied a lower trophic level, there was a 7.6% decrease in community-level native responses over the typical range in invasive abundance investigated in the literature

##########################################################

## Figure 3a: Population response for invader at higher trophic #######################################


par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res4a[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Population response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res4a[,1],res4a[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res4a[,3],type="l", lwd = 2)


## Figure 3b: Population response for invader at same trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res4b[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Population response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res4b[,1],res4b[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res4b[,3],type="l", lwd = 2)



## Figure 3c: Population response for invader at lower trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res4c[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Population response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res4c[,1],res4c[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res4c[,3],type="l", lwd = 2)



## Figure 3d: Community response for invader at higher trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res4d[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Community response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res4d[,1],res4d[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res4d[,3],type="l", lwd = 2)


## Figure 3e: Community response for invader at same trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res4e[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Community response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res4e[,1],res4e[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res4e[,3],type="l", lwd = 2)



## Figure 3f: Community response for invader at lower trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res4f[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Community response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res4f[,1],res4f[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res4f[,3],type="l", lwd = 2)



##################################################################
## Figure 4: Slopes by diversity metric 
##################################################################
avi.ES4.s <- avi.ES.s[which(avi.ES.s$Multi_spp_resp == "MULTIPLE"),]
avi.ES4.s <- avi.ES4.s[-c(which(avi.ES4.s$Response_type=="Abundance"),which(avi.ES4.s$Response_type=="Other")),]
# remove empty levels (clean up)
avi.ES4.s <- droplevels (avi.ES4.s)


ms29a <- MCMCglmm(Int.u ~ Response_type, random = ~ Article_ID, mev = avi.ES4.s$Int.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES4.s)
ms29b <- MCMCglmm(Lin.u ~ Response_type, random = ~ Article_ID, mev = avi.ES4.s$Lin.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES4.s)
ms29c <- MCMCglmm(Pol.u ~ Response_type, random = ~ Article_ID, mev = avi.ES4.s$Pol.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES4.s)

summary(ms29b)
# The linear term for diversity responses is significantly less than zero p < 0.001
# The linear term for evenness is not significantly different from that for diversity  p = 0.278
# The linear term for richness is significantly greater than that for diversity  p = 0.036

# set reference level of Trophic_level to intra, to examine intra vs. below
avi.ES4.s$Response_type <- relevel (avi.ES4.s$Response_type, ref = 3)

#model m29b2
m29b2 <- MCMCglmm(Lin.u ~ Response_type, random = ~ Article_ID, mev = avi.ES4.s$Lin.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES4.s)
summary(m29b2)

# Linear effect size for evenness is significantly lower than (more negative) than that for richness below p = 0.004

# VERY IMPORTANT: Reset the Trophic_level variable to original indexing
avi.ES3.s$Trophic_level <- relevel (avi.ES3.s$Trophic_level, ref = 2)


summary(ms29c)
# The polynomial term for diversity responses is not significantly different from zero p = 0.888
# The polynomial term for evenness is not significantly different from that for diversity  p = 0.200
# The polynomial term for richness is not significantly different from that for diversity  p = 0.188


# set reference level of Trophic_level to intra, to examine intra vs. below
avi.ES4.s$Response_type <- relevel (avi.ES4.s$Response_type, ref = 2)

#model m29c2
m29c2 <- MCMCglmm(Pol.u ~ Response_type, random = ~ Article_ID, mev = avi.ES4.s$Pol.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES4.s)
summary(m29c2)

# Polynomial effect size for evenness is significantly greater than that for richness p = 0.012

# VERY IMPORTANT: Reset the Trophic_level variable to original indexing
avi.ES3.s$Trophic_level <- relevel (avi.ES3.s$Trophic_level, ref = 2)




ms29.sol <- cbind(ms29a$Sol[, "(Intercept)"] + ms29a$Sol[, "Response_typeRichness"],
				ms29b$Sol[, "(Intercept)"] + ms29b$Sol[, "Response_typeRichness"],
				ms29c$Sol[, "(Intercept)"] + ms29c$Sol[, "Response_typeRichness"],
				ms29a$Sol[, "(Intercept)"], 
				ms29b$Sol[, "(Intercept)"], 
				ms29c$Sol[, "(Intercept)"], 
				ms29a$Sol[, "(Intercept)"]+ ms29a$Sol[, "Response_typeEvenness"],
				ms29b$Sol[, "(Intercept)"]+ ms29b$Sol[, "Response_typeEvenness"],
				ms29c$Sol[, "(Intercept)"]+ ms29c$Sol[, "Response_typeEvenness"])


# Estimates for 9 regression term effect sizes
apply(mcmc(ms29.sol), 2, mean)
# Ordering of estimates is as follows, these contain information to add to plots
# 1 Richness, Intercept
# 2 Richness, Linear
# 3 Richness, Polynomial
# 4 Diversity, Intercept
# 5 Diversity, Linear
# 6 Diversity, Polynomial
# 7 Evenness, Intercept
# 8 Evenness, Linear
# 9 Evenness, Polynomial


# p-values for 9 regression term effect sizes, comparing to zero
apply(mcmc(ms29.sol), 2, pMCMC)
# Ordering is as in comments immediately preceding, some of these p-values to be added to plots as star symbols



# split by row, for lapply convenience later...
ms29.sol <- split(ms29.sol,1:1000)
ix <- (-500:500)/1000

####### (i) Richness: predictions and credible zone #############################

# matrix to catch results for 95% credible zone
res29a <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms29.sol, function(y){
    ny <- y[3]*(ix[i]^2) + y[2]*ix[i] + y[1]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res29a[i,] <- c(pred1.int,pred1.u)
}

####### (ii) Diversity: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res29b <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms29.sol, function(y){
	# Note indexing within the square brackets refers to correct columns in ms29.sol
    ny <- y[6]*(ix[i]^2) + y[5]*ix[i] + y[4]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res29b[i,] <- c(pred1.int,pred1.u)
}

####### (iii) Evenness, Below: predictions and credible zone ##########################

# matrix to catch results for 95% credible zone
res29c <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms29.sol, function(y){
	# Note indexing within the square brackets refers to correct columns in ms29.sol
    ny <- y[9]*(ix[i]^2) + y[8]*ix[i] + y[7]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res29c[i,] <- c(pred1.int,pred1.u)
}

######## Proportional change in native responses ###################################
######## over the typical range in invasive species' ###############################
######## abundance investigated in the literature ##################################
# note: changes are on average over meta-analysed studies, the response in individual studies varies.

res29a[1000,3] - res29a[1,3]
# There was a 10.9% decrease in native species richness over the typical range in invasive abundance investigated in the literature
res29b[1000,3] - res29b[1,3]
# There was a 23.4% decrease in native species diversity over the typical range in invasive abundance investigated in the literature
res29c[1000,3] - res29c[1,3]
# There was a 29.8% decrease in native species evenness over the typical range in invasive abundance investigated in the literature

##########################################################

## Figure 4a: Slopes for richness analyses #######################################


par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
par(pty="s") #square!!
plot((ix+0.5),res29a[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res29a[,1],res29a[length(ix):1,2]),border = NA, col = "red")
lines((ix+0.5),res29a[,3],type="l", lwd = 2)



###

## Figure 4b: Slopes for diversity analyses #######################################


par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
par(pty="s") #square!!
plot((ix+0.5),res29b[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2, new = F)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res29b[,1],res29b[length(ix):1,2]),border = NA, col = "cyan")
lines((ix+0.5),res29b[,3],type="l", lwd = 2)




## Figure 4c: Slopes for evenness analyses #######################################


par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
par(pty="s") #square!!
plot((ix+0.5),res29c[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res29c[,1],res29c[length(ix):1,2]),border = NA, col = "blue")
lines((ix+0.5),res29c[,3],type="l", lwd = 2)


##################################################################
## Analysis S3.3: Recipient habitat/ ecosystem by trophic level 
##################################################################



ms10a <- MCMCglmm(Int.u ~ Inv_habitat + Trophic_level, random = ~ Article_ID, mev = avi.ES3.s$Int.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3.s)
ms10b <- MCMCglmm(Lin.u ~ Inv_habitat + Trophic_level, random = ~ Article_ID, mev = avi.ES3.s$Lin.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3.s)
ms10c <- MCMCglmm(Pol.u ~ Inv_habitat + Trophic_level, random = ~ Article_ID, mev = avi.ES3.s$Pol.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3.s)

summary(ms10b)
# The linear term for Aquatic studies is significantly less than zero p < 0.001
# Neither the linear term for terrestrial nor marine studies are significantly different from that for aquatic studies (p = 0.328; p = 0.376)

# set reference level of Inv_habitat to terrestrial, to examine terr vs. marine
avi.ES3.s$Inv_habitat <- relevel (avi.ES3.s$Inv_habitat, ref = 3)

#model m10b2
m10b2 <- MCMCglmm(Lin.u ~ Inv_habitat + Trophic_level, random = ~ Article_ID, mev = avi.ES3.s$Lin.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3.s)
summary(m10b2)

# Linear effect size for terrestrial studies is significantly greater than (less negative) than that for  marine studies p = 0.024

# VERY IMPORTANT: Reset the Inv_habitat variable to original indexing
avi.ES3.s$Inv_habitat <- relevel (avi.ES3.s$Inv_habitat, ref = 2)

summary(ms10c)
# The polynomial term for Aquatic studies is significantly less than zero p < 0.001
# Neither the polynomial term for terrestrial nor marine studies are significantly different from that for aquatic studies (p = 0.112; p = 0.152)

# set reference level of Inv_habitat to terrestrial, to examine terr vs. marine
avi.ES3.s$Inv_habitat <- relevel (avi.ES3.s$Inv_habitat, ref = 3)

#model m10c2
m10c2 <- MCMCglmm(Pol.u ~ Inv_habitat + Trophic_level, random = ~ Article_ID, mev = avi.ES3.s$Pol.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3.s)
summary(m10c2)

# Polynomial effect size for terrestrial studies not significantly different from that for  marine studies p = 0.878


# VERY IMPORTANT: Reset the Inv_habitat variable to original indexing
avi.ES3.s$Inv_habitat <- relevel (avi.ES3.s$Inv_habitat, ref = 2)


ms10.sol <- cbind(ms10a$Sol[, "(Intercept)"] + ms10a$Sol[, "Inv_habitatterrestrial"],
				ms10b$Sol[, "(Intercept)"] + ms10b$Sol[, "Inv_habitatterrestrial"],
				ms10c$Sol[, "(Intercept)"] + ms10c$Sol[, "Inv_habitatterrestrial"],

				ms10a$Sol[, "(Intercept)"], 
				ms10b$Sol[, "(Intercept)"], 
				ms10c$Sol[, "(Intercept)"], 

				ms10a$Sol[, "(Intercept)"]+ ms10a$Sol[, "Inv_habitatmarine"],
				ms10b$Sol[, "(Intercept)"]+ ms10b$Sol[, "Inv_habitatmarine"],
				ms10c$Sol[, "(Intercept)"]+ ms10c$Sol[, "Inv_habitatmarine"],

				ms10a$Sol[, "(Intercept)"] + ms10a$Sol[, "Inv_habitatterrestrial"] + ms10a$Sol[, "Trophic_levelIntra"], 
				ms10b$Sol[, "(Intercept)"] + ms10b$Sol[, "Inv_habitatterrestrial"] + ms10b$Sol[, "Trophic_levelIntra"], 
				ms10c$Sol[, "(Intercept)"] + ms10c$Sol[, "Inv_habitatterrestrial"] + ms10c$Sol[, "Trophic_levelIntra"], 

				ms10a$Sol[, "(Intercept)"] + ms10a$Sol[, "Trophic_levelIntra"], 
				ms10b$Sol[, "(Intercept)"] + ms10b$Sol[, "Trophic_levelIntra"], 
				ms10c$Sol[, "(Intercept)"] + ms10c$Sol[, "Trophic_levelIntra"], 

				ms10a$Sol[, "(Intercept)"]+ ms10a$Sol[, "Inv_habitatmarine"] + ms10a$Sol[, "Trophic_levelIntra"], 
				ms10b$Sol[, "(Intercept)"]+ ms10b$Sol[, "Inv_habitatmarine"] + ms10b$Sol[, "Trophic_levelIntra"], 
				ms10c$Sol[, "(Intercept)"]+ ms10c$Sol[, "Inv_habitatmarine"] + ms10c$Sol[, "Trophic_levelIntra"], 

				ms10a$Sol[, "(Intercept)"] + ms10a$Sol[, "Inv_habitatterrestrial"] + ms10a$Sol[, "Trophic_levelBelow"], 
				ms10b$Sol[, "(Intercept)"] + ms10b$Sol[, "Inv_habitatterrestrial"] + ms10b$Sol[, "Trophic_levelBelow"], 
				ms10c$Sol[, "(Intercept)"] + ms10c$Sol[, "Inv_habitatterrestrial"] + ms10c$Sol[, "Trophic_levelBelow"], 

				ms10a$Sol[, "(Intercept)"] + ms10a$Sol[, "Trophic_levelBelow"], 
				ms10b$Sol[, "(Intercept)"] + ms10b$Sol[, "Trophic_levelBelow"], 
				ms10c$Sol[, "(Intercept)"] + ms10c$Sol[, "Trophic_levelBelow"], 

				ms10a$Sol[, "(Intercept)"]+ ms10a$Sol[, "Inv_habitatmarine"] + ms10a$Sol[, "Trophic_levelBelow"], 
				ms10b$Sol[, "(Intercept)"]+ ms10b$Sol[, "Inv_habitatmarine"] + ms10b$Sol[, "Trophic_levelBelow"], 
				ms10c$Sol[, "(Intercept)"]+ ms10c$Sol[, "Inv_habitatmarine"] + ms10c$Sol[, "Trophic_levelBelow"]) 


# Estimates for 27 regression term effect sizes
apply(mcmc(ms10.sol), 2, mean)
# Ordering of estimates is as follows, these contain information to add to plots
# 1 Terrestrial, Above, Intercept
# 2 Terrestrial, Above, Linear
# 3 Terrestrial, Above, Polynomial
# 4 Aquatic, Above, Intercept
# 5 Aquatic, Above, Linear
# 6 Aquatic, Above, Polynomial
# 7 Marine, Above, Intercept
# 8 Marine, Above, Linear
# 9 Marine, Above, Polynomial
# 10 Terrestrial, Intra, Intercept
# 11 Terrestrial, Intra, Linear
# 12 Terrestrial, Intra, Polynomial
# 13 Aquatic, Intra, Intercept
# 14 Aquatic, Intra, Linear
# 15 Aquatic, Intra, Polynomial
# 16 Marine, Intra, Intercept
# 17 Marine, Intra, Linear
# 18 Marine, Intra, Polynomial
# 19 Terrestrial, Below, Intercept
# 20 Terrestrial, Below, Linear
# 21 Terrestrial, Below, Polynomial
# 22 Aquatic, Below, Intercept
# 23 Aquatic, Below, Linear
# 24 Aquatic, Below, Polynomial
# 25 Marine, Below, Intercept
# 26 Marine, Below, Linear
# 27 Marine, Below, Polynomial


# p-values for 27 regression term effect sizes, comparing to zero
apply(mcmc(ms10.sol), 2, pMCMC)
# Ordering is as in comments immediately preceding, some of these p-values to be added to plots as star symbols



# split by row, for lapply convenience later...
ms10.sol <- split(ms10.sol,1:1000)
ix <- (-500:500)/1000

####### (i) Terrestrial, above: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res10a <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms10.sol, function(y){
    ny <- y[3]*(ix[i]^2) + y[2]*ix[i] + y[1]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res10a[i,] <- c(pred1.int,pred1.u)
}

####### (ii) Aquatic, above: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res10b <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms10.sol, function(y){
    ny <- y[6]*(ix[i]^2) + y[5]*ix[i] + y[4]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res10b[i,] <- c(pred1.int,pred1.u)
}

####### (iii) Marine, above: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res10c <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms10.sol, function(y){
    ny <- y[9]*(ix[i]^2) + y[8]*ix[i] + y[7]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res10c[i,] <- c(pred1.int,pred1.u)
}

####### (iv) Terrestrial, intra: predictions and credible zone #########################

# matrix to catch results for 95% credible zone
res10d <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms10.sol, function(y){
    ny <- y[12]*(ix[i]^2) + y[11]*ix[i] + y[10]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res10d[i,] <- c(pred1.int,pred1.u)
}

####### (v) Aquatic, intra: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res10e <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms10.sol, function(y){
    ny <- y[15]*(ix[i]^2) + y[14]*ix[i] + y[13]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res10e[i,] <- c(pred1.int,pred1.u)
}

####### (vi) Marine, intra: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res10f <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms10.sol, function(y){
    ny <- y[18]*(ix[i]^2) + y[17]*ix[i] + y[16]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res10f[i,] <- c(pred1.int,pred1.u)
}

####### (vii) Terrestrial, below: predictions and credible zone #########################

# matrix to catch results for 95% credible zone
res10g <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms10.sol, function(y){
    ny <- y[21]*(ix[i]^2) + y[20]*ix[i] + y[19]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res10g[i,] <- c(pred1.int,pred1.u)
}

####### (viii) Aquatic, below: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res10h <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms10.sol, function(y){
    ny <- y[24]*(ix[i]^2) + y[23]*ix[i] + y[22]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res10h[i,] <- c(pred1.int,pred1.u)
}

####### (ix) Marine, below: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res10i <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms10.sol, function(y){
    ny <- y[27]*(ix[i]^2) + y[26]*ix[i] + y[25]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res10i[i,] <- c(pred1.int,pred1.u)
}

######## Proportional change in native responses ###################################
######## over the typical range in invasive species' ###############################
######## abundance investigated in the literature ##################################
# note: changes are on average over meta-analysed studies, the response in individual studies varies.

res10a[1000,3] - res10a[1,3]
# Where the invasive species was at a higher trophic level, there was a 54.9% decrease in terrestrial native responses over the typical range in invasive abundance investigated in the literature
res10b[1000,3] - res10b[1,3]
# Where the invasive species was at a higher trophic level, there was a 45.7% decrease in aquatic native responses over the typical range in invasive abundance investigated in the literature
res10c[1000,3] - res10c[1,3]
# Where the invasive species was at a higher trophic level, there was a 35.0% decrease in marine native responses over the typical range in invasive abundance investigated in the literature
res10d[1000,3] - res10d[1,3]
# Where the invasive species was at the same trophic level, there was a 30.0% decrease in terrestrial native responses over the typical range in invasive abundance investigated in the literature
res10e[1000,3] - res10e[1,3]
# Where the invasive species was at the same trophic level, there was a 20.8% decrease in aquatic native responses over the typical range in invasive abundance investigated in the literature
res10f[1000,3] - res10f[1,3]
# Where the invasive species was at the same trophic level, there was a 10.1% decrease in marine native responses over the typical range in invasive abundance investigated in the literature
res10g[1000,3] - res10g[1,3]
# Where the invasive species was at a lower trophic level, there was a 11.4% decrease in terrestrial native responses over the typical range in invasive abundance investigated in the literature
res10h[1000,3] - res10h[1,3]
# Where the invasive species was at a lower trophic level, there was a 2.2% decrease in aquatic native responses over the typical range in invasive abundance investigated in the literature
res10i[1000,3] - res10i[1,3]
# Where the invasive species was at a lower trophic level, there was a 8.5% increase in marine native responses over the typical range in invasive abundance investigated in the literature


##########################################################




## Figure S3.3a: Terrestrial response for invader at higher trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res10a[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response (Terr)", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res10a[,1],res10a[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res10a[,3],type="l", lwd = 2)



## Figure S3.3d: Aquatic response for invader at higher trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res10b[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response (Aqua)", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res10b[,1],res10b[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res10b[,3],type="l", lwd = 2)



## Figure S3.3g: Marine response for invader at higher trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res10c[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response (Mar)", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res10c[,1],res10c[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res10c[,3],type="l", lwd = 2)



## Figure S3.3b: Terrestrial response for invader at same trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res10d[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response (Terr)", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res10d[,1],res10d[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res10d[,3],type="l", lwd = 2)


## Figure S3.3e: Aquatic response for invader at same trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res10e[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response (Aqua)", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res10e[,1],res10e[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res10e[,3],type="l", lwd = 2)


## Figure S3.3h: Marine response for invader at same trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res10f[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response (Mar)", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res10f[,1],res10f[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res10f[,3],type="l", lwd = 2)


## Figure S3.3c: Terrestrial response for invader at lower trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res10g[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response (Terr)", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res10g[,1],res10g[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res10g[,3],type="l", lwd = 2)


## Figure S3.3f: Aquatic response for invader at lower trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res10h[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response (Aqua)", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res10h[,1],res10h[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res10h[,3],type="l", lwd = 2)


## Figure S3.3i: Marine response for invader at lower trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res10i[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res10i[,1],res10i[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res10i[,3],type="l", lwd = 2)




##################################################################
## Analysis S3.4: Invasive plants vs. animals by trophic level 
##################################################################

avi.ES5.s <- avi.ES3.s[-which(avi.ES3.s$Inv_kingdom =="Bacteria"),]
avi.ES5.s <- droplevels(avi.ES5.s)


ms19a <- MCMCglmm(Int.u ~ Inv_kingdom + Trophic_level, random = ~ Article_ID, mev = avi.ES5.s$Int.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES5.s)
ms19b <- MCMCglmm(Lin.u ~ Inv_kingdom + Trophic_level, random = ~ Article_ID, mev = avi.ES5.s$Lin.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES5.s)
ms19c <- MCMCglmm(Pol.u ~ Inv_kingdom + Trophic_level, random = ~ Article_ID, mev = avi.ES5.s$Pol.s, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES5.s)

summary(ms19b) 
# Linear effect size for invasive animals is significantly less than zero (p < 0.001)
# Linear effect sizes of invasive plants and animals do not differ (p = 0.186)

summary(ms19c) 
# Polynomial effect size for invasive animals is significantly greater than zero (p < 0.001)
# Polynomial effect sizes significantly lower in invasive plants than invasive animals (p = 0.036)


ms19.sol <- cbind(ms19a$Sol[, "(Intercept)"],
				ms19b$Sol[, "(Intercept)"],
				ms19c$Sol[, "(Intercept)"],
				ms19a$Sol[, "(Intercept)"] + ms19a$Sol[, "Trophic_levelIntra"],
				ms19b$Sol[, "(Intercept)"] + ms19b$Sol[, "Trophic_levelIntra"],
				ms19c$Sol[, "(Intercept)"] + ms19c$Sol[, "Trophic_levelIntra"],
				ms19a$Sol[, "(Intercept)"] + ms19a$Sol[, "Trophic_levelBelow"],
				ms19b$Sol[, "(Intercept)"] + ms19b$Sol[, "Trophic_levelBelow"],
				ms19c$Sol[, "(Intercept)"] + ms19c$Sol[, "Trophic_levelBelow"],
				ms19a$Sol[, "(Intercept)"] + ms19a$Sol[, "Inv_kingdomPlant"],
				ms19b$Sol[, "(Intercept)"] + ms19b$Sol[, "Inv_kingdomPlant"],
				ms19c$Sol[, "(Intercept)"] + ms19c$Sol[, "Inv_kingdomPlant"],
				ms19a$Sol[, "(Intercept)"] + ms19a$Sol[, "Trophic_levelIntra"]+ ms19a$Sol[, "Inv_kingdomPlant"],
				ms19b$Sol[, "(Intercept)"] + ms19b$Sol[, "Trophic_levelIntra"]+ ms19b$Sol[, "Inv_kingdomPlant"],
				ms19c$Sol[, "(Intercept)"] + ms19c$Sol[, "Trophic_levelIntra"]+ ms19c$Sol[, "Inv_kingdomPlant"],
				ms19a$Sol[, "(Intercept)"] + ms19a$Sol[, "Trophic_levelBelow"]+ ms19a$Sol[, "Inv_kingdomPlant"],
				ms19b$Sol[, "(Intercept)"] + ms19b$Sol[, "Trophic_levelBelow"]+ ms19b$Sol[, "Inv_kingdomPlant"],
				ms19c$Sol[, "(Intercept)"] + ms19c$Sol[, "Trophic_levelBelow"]+ ms19c$Sol[, "Inv_kingdomPlant"])


# Estimates for 18 regression term effect sizes
apply(mcmc(ms19.sol), 2, mean)
# Ordering of estimates is as follows, these contain information to add to plots
# 1 Animal, Above, Intercept
# 2 Animal, Above, Linear
# 3 Animal, Above, Polynomial
# 4 Animal, Intra, Intercept
# 5 Animal, Intra, Linear
# 6 Animal, Intra, Polynomial
# 7 Animal, Below, Intercept
# 8 Animal, Below, Linear
# 9 Animal, Below, Polynomial
# 10 Plant, Above, Intercept
# 11 Plant, Above, Linear
# 12 Plant, Above, Polynomial
# 13 Plant, Intra, Intercept
# 14 Plant, Intra, Linear
# 15 Plant, Intra, Polynomial
# 16 Plant, Below, Intercept
# 17 Plant, Below, Linear
# 18 Plant, Below, Polynomial



# p-values for 18 regression term effect sizes, comparing to zero
apply(mcmc(ms19.sol), 2, pMCMC)
# Ordering is as in comments immediately preceding, some of these p-values to be added to plots as star symbols



# split by row, for lapply convenience later...
ms19.sol <- split(ms19.sol,1:1000)
ix <- (-500:500)/1000

####### (i) Animal, above: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res19a <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms19.sol, function(y){
    ny <- y[3]*(ix[i]^2) + y[2]*ix[i] + y[1]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res19a[i,] <- c(pred1.int,pred1.u)
}

####### (ii) Animal, intra: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res19b <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms19.sol, function(y){
    ny <- y[6]*(ix[i]^2) + y[5]*ix[i] + y[4]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res19b[i,] <- c(pred1.int,pred1.u)
}

####### (iii) Animal, below: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res19c <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms19.sol, function(y){
    ny <- y[9]*(ix[i]^2) + y[8]*ix[i] + y[7]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res19c[i,] <- c(pred1.int,pred1.u)
}

####### (iv) Plant, above: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res19d <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms19.sol, function(y){
    ny <- y[12]*(ix[i]^2) + y[11]*ix[i] + y[10]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res19d[i,] <- c(pred1.int,pred1.u)
}

####### (v) Plant, intra: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res19e <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms19.sol, function(y){
    ny <- y[15]*(ix[i]^2) + y[14]*ix[i] + y[13]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res19e[i,] <- c(pred1.int,pred1.u)
}

####### (vi) Plant, below: predictions and credible zone ###########################

# matrix to catch results for 95% credible zone
res19f <- matrix(NA,length(ix),3)

for (i in 1:length(ix)){
  pred1 <- lapply(ms19.sol, function(y){
    ny <- y[18]*(ix[i]^2) + y[17]*ix[i] + y[16]
    return (ny)
  })
  pred1.int <- HPDinterval(mcmc(unlist(pred1)))
  
  #	expected functional curve is summarised using the mean of predicted posterior samples
  pred1.u <- mean(unlist(pred1))
  res19f[i,] <- c(pred1.int,pred1.u)
}


######## Proportional change in native responses ###################################
######## over the typical range in invasive species' ###############################
######## abundance investigated in the literature ##################################
# note: changes are on average over meta-analysed studies, the response in individual studies varies.

res19a[1000,3] - res19a[1,3]
# Where the invasive species was at a higher trophic level, there was a 47.3% decrease in  native animal responses over the typical range in invasive abundance investigated in the literature
res19b[1000,3] - res19b[1,3]
# Where the invasive species was at the same trophic level, there was a 18.2% decrease in native animal responses over the typical range in invasive abundance investigated in the literature
res19c[1000,3] - res19c[1,3]
# Where the invasive species was at a lower trophic level, there was a 1.3% increase in native animal responses over the typical range in invasive abundance investigated in the literature
res19d[1000,3] - res19d[1,3]
# Where the invasive species was at a higher trophic level, there was a 58.0% decrease in native plant responses over the typical range in invasive abundance investigated in the literature
res19e[1000,3] - res19e[1,3]
# Where the invasive species was at the same trophic level, there was a 28.9% decrease in native plant responses over the typical range in invasive abundance investigated in the literature
res19f[1000,3] - res19f[1,3]
# Where the invasive species was at a lower trophic level, there was a 9.4% decrease in native plant responses over the typical range in invasive abundance investigated in the literature

##########################################################



## Figure S3.4a: Invasive animal at higher trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res19a[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res19a[,1],res19a[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res19a[,3],type="l", lwd = 2)



## Figure S3.4b: Invasive animal at same trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res19b[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res19b[,1],res19b[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res19b[,3],type="l", lwd = 2)


## Figure S3.4c: Invasive animal at lower trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res19c[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res19c[,1],res19c[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res19c[,3],type="l", lwd = 2)


## Figure S3.4d: Invasive plant at same trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res19e[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res19e[,1],res19e[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res19e[,3],type="l", lwd = 2)


## Figure S3.4e: Invasive plant at lower trophic #######################################

par(mfrow=c(1,1))
par(mar = c(5,5,2,2))
plot((ix+0.5),res19f[,3],type="n",xlim=c(0,1),ylim=c(0,1), ylab = "Native response", xlab = "Invasive abundance", cex.lab = 2, cex.axis = 2)
polygon(x = c(ix,ix[length(ix):1])+0.5, y = c(res19f[,1],res19f[length(ix):1,2]),border = NA, col = "darkgrey")
lines((ix+0.5),res19f[,3],type="l", lwd = 2)





##################################################################
##################################################################
##################################################################
# SECTION 6. function to look see how the sign of the raw data is distributed, both for invasive and native, run the hashed code following the function


# loc.xy <- function (data) {
# 	
# 	xx <- split(data,list(data$Article_ID,data$Study_ID),drop=T)	
# yy <- lapply(xx, function(y){
# 	c(sign(min(y$Response)), sign(max(y$Response)), sign(min(y$Abundance_Invader)), sign(max(y$Abundance_Invader)))
# })
# 
# yy2 <- do.call(rbind,yy)
# 
# return (yy2)
# 	
# }
# 
# test <- loc.xy(data)
# table(test[,1],test[,2]) #min, max signs for native response
# table(test[,3],test[,4]) #min, max signs for invader abundance
##################################################################
##################################################################
##################################################################

