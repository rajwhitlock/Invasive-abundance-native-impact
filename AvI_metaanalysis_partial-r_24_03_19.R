#Meta-analysis code and analysis framework written by Raj Whitlock (RW) 05/18/2017
#Code last updated by Bethany Bradley & RW 03/24/2019
#Overview: this code allows raw data pre-processing to calculate "partial-r" effect sizes, meta-analysis of processed effect-size data and figure creation from model output
#You will need raw input files "AvI_data_sheet.txt" and "AvI_attribute_sheet.txt" as well as the R script 'Bar_plot_error_bars.R'. Place the datasets in your working directory
# n.b. all forest plots are bivariate forest plots summarising both linear and polynomial components to invasive abundance – native response relationships

##################################################################
##################################################################
source("/Directory/path/to/script/Bar_plot_error_bars.R")

setwd("/Directory/path/to/working_directory")

#load necessary libraries
library(MASS)
library(MCMCglmm)
library(metafor)
library(dplyr)

# Code to create effect size dataset from raw data

##################################################################
##################################################################
##################################################################

# SECTION 1. The data
full_data <- read.delim("AvI_data_sheet.txt")
#Read in the AvI data sheet containing all Article and Study level attributes associated with each observation
aie_db <- read.delim("AvI_attribute_sheet.txt", header=T)

##################################################################
##################################################################
##################################################################

# SECTION 2. Needed utility functions:
# Equation from Nakagawa & Cuthill 2007 in Biol Reviews 82:591.  
# This function calculates a partial correlation coefficient (partial.r) for effect size generation, later in the script, from a t-value, and df

partial.r<-function(t.val,df){
  r<-t.val/sqrt((t.val)^2+df)
  names(r)<-"effect size r"
  return(r)}

# This is a function to calculate approximate Bayesian p-vales from posterior samples from MCMCglmm models, following Hadfield JD (2010) MCMC methods for multi-response generalized linear mixed models: the MCMCglmm R package. Journal of Statistical Software 33, 1–22.

pMCMC <- function(data){
res <- 2*(min(length(which(data >0)),length(which(data <0))))/length(data)
if (res<0.001){res <- "<0.001"}
return(as.character(res)) 
}

##################################################################
##################################################################
##################################################################

# SECTION 3. Data pre-processing to calculate effect sizes


#Define some "columns to pass" (c2p) through data pre-processing that are cut out and joined back into the effect size dataset prior to analysis
c2p <- c(which(names(full_data)=="Article_ID"),which(names(full_data)=="Study_ID"),which(names(full_data)=="Art_stud"))



##################################################################
##################################################################
####################START OF FUNCTION CREATING "PARTIAL-R" EFFECT SIZES########

# The following function will process the raw data. The argument columns.to.pass specifies article-level variates that should be retained, excluding the unique article identifier, which will be retained automatically

avi.es.calc <- function (data,columns.to.pass, rsc.y = T, rsc.x = T){ 
  
  avi.cov <- data[,c2p]
  avi.cov2 <- split(avi.cov,f = list(avi.cov[,1],avi.cov[,2]),drop=T)
  avi.cov3 <- lapply(avi.cov2,function(x){
    x[1,]
  })
  avi.cov <- do.call(rbind,avi.cov3)
  
  xx <- split(data,list(data$Article_ID,data$Study_ID),drop=T)	
  yy <- lapply(xx, function(y){
    
    #Centering each invader abundance by StudyID is necessary for analysis/interpretability
    #Schielzeth 2010 in Methods Ecol & Evol
    #Centering reduces the dependency between slope and curvature estimates
    
    n <- dim(y)[1]
    
    if (rsc.x == T){
      if(sign (min (y$Abundance_Invader))==-1){offset.x <- min (y$Abundance_Invader)} else {offset.x <- 0}
      inv.recentred <- (y$Abundance_Invader - offset.x)/(max(y$Abundance_Invader) - offset.x)
      inv.recentred <- inv.recentred - mean(inv.recentred)
    } else {
      inv.recentred <- y$Abundance_Invader - mean(y$Abundance_Invader)
    }
    
    # rescaling the data between 0-1 (if rsc.y is True)	
    if (rsc.y == T){
      if(sign (min (y$Response))==-1){offset.y <- min (y$Response)} else {offset.y <- 0}
      y.resp <- (y$Response - offset.y) / (max(y$Response) - offset.y)
    } else {y.resp <- y$Response}
    
    #mev is assessing measurement error variance based on the number of rows in each StudyID
    #more rows, less variance
    mev <- 1/(n-3)
    
    #If else statement making sure the studies have >3 unique sets of values
    if (n > 3 && length(unique(y$Abundance_Invader)) > 3){
      
      #m1, below is a linear model of the response, which includes both the linear and polynomial components of the model fit    
      
      m1 <- lm(y$Response ~ poly(inv.recentred,2, raw = T)) #model in PNAS paper
      	# Raw polynomials were chosen so that a consistently scaled model matrix applied in each of the studies, allowing comparability and synthesis of effect sizes from different studies
      a.1 <- summary(m1)$coef[2,3] #t value of linear term
      a.2 <- summary(m1)$coef[3,3] #t value of polynomial term
      
      r.1 <- partial.r(a.1, df = 1)		# create partial correlation coef using t value corresponding to coefficient for linear term
      r.2 <- partial.r(a.2, df = 1)   # create partial correlation coef using t value corresponding to coefficient for polynomial term
      c.m1 <- cov2cor(vcov(m1))[3,2]  # calculate within-study correlation of linear and polynomial predictors
      
      
      # Effect size is fisher z transformation of partial correlation coefficient:
      es.1 <- 0.5*log((1 + r.1)/(1 - r.1))
      es.2 <- 0.5*log((1 + r.2)/(1 - r.2))
      
      #Repeating analysis for studies with 3 abundance categories with only the linear fit (no polynomial)
    }   else    {
      

      m1 <- lm(y.resp ~ inv.recentred)
      a.1 <- summary(m1)$coef[2,3]
      
      r.1 <- partial.r(a.1, df = 1)		# create partial correlation coef using t value corresponding to coefficient for quadratic term
      r.2 <- NA
      c.m1 <- NA
      # Effect size is fisher z transformation of correlation coefficient:
      
      es.1 <- 0.5*log((1 + r.1)/(1 - r.1))
      es.2 <- NA
    }
    
    return(data.frame("Article_ID" = as.character(y$Article_ID[1]), "Study_ID" = as.character(y$Study_ID[1]), "ES.linear" = es.1, "ES.poly" = es.2, "mev" = mev, "cor1.2" = c.m1))
  })
  
  #binding together all of the individual rows into a master sheet of article, study & linear/poly estimates
  yy2 <- do.call(rbind,yy)
  
  #binding the year column back in at the end
  res <- merge (avi.cov,yy2)
  return (res)
}

##################################################################
####################END OF FUNCTION CREATING PARTIAL.R EFFECT SIZES########
##################################################################

# SECTION 4. Now run the function and create the effect size dataset:

avi.ES <- avi.es.calc(full_data,columns.to.pass=c2p, rsc.x = T, rsc.y = T) 

#database internal join of effect estimates with article level covariates
avi.ES <- merge(avi.ES, aie_db)

#Four studies removed: BRAN2010_1, BRAN2010_3 and TERA2011_2 all have <3 response values so poly term can't be calculated
#TERA2011_4 poly.u = -464 (only five lines in study) - large outlier poly value.  
avi.ES <- avi.ES[-which(avi.ES$Art_stud=="BRAN2010_1"),]
avi.ES <- avi.ES[-which(avi.ES$Art_stud=="TERA2011_2"),]
avi.ES <- avi.ES[-which(avi.ES$Art_stud=="BRAN2010_3"),]
avi.ES <- avi.ES[-which(avi.ES$Art_stud=="TERA2011_4"),]

#dim(avi.ES)	# Should be 1258 rows, 31 columns

##################################################################
##################################################################

# SECTION 5. Meta-analysis. The code in this section specifies models that can be run on the effect size data, and code for producing figures from model output. Take care, the correspondence between figure numbers indicated in this code and figure numbers in the manuscript and its supplementary materials may need checking

# We will use a Bayesian random effects model with article/ paper level random effects

#uniform prior on the standard deviation of the random effects, for both residual variance (R structure), and study level variance (G structure)

prior <- list(R=list(V = 1e-10, nu = -1), G = list(G1 = list(V = 1e-10, nu = -1)))

##################################################################
##################################################################
## Analysis 1: Global meta-analysis ## (***not included in paper figures***)

# random effects meta-analysis with additional random effects for study, no fixed effects bar the intercept
# RW adapted the burnin and thin numbers based on the model summary and sample sizes
m1a <- MCMCglmm(ES.linear ~ 1, random = ~ Article_ID, mev = avi.ES$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES)
#linear model direction summary for all StudyIDs
summary(m1a)

m1b <- MCMCglmm(ES.poly ~ 1, random = ~ Article_ID, mev = avi.ES$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES)
summary(m1b)

#################### Bivariate Forest Plot Figure Global ###############################
res1a <- summary(m1a)$solutions #linear
res1b <- summary(m1b)$solutions #poly

res1a_1b <- rbind (res1b,res1a)
#sample size:
n1 <- length(unique(avi.ES$Article_ID))		# n = 202 articles
n2 <- length(unique(avi.ES$Art_stud))		# n = 1258 studies

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
plot(res1a_1b[,1],c(0.5:1.5),xlim=c(-2,1.5),pch=22,bg=c(rep("black",2)),yaxt="n",bty="n",ylab="",xlab="Effect size",
     cex.lab=1.3,cex.axis=1.3,xaxp=c(-1,0.5,3),ylim=c(0,2))
axis(side=3,labels=F,at=c(-1,-0.5,-0,0.5))
e.bars.points2(res1a_1b[,1], z1 = res1a_1b[,3], z2 = res1a_1b[,2], x = c(0.5,1.5),cap=0.1,clr=c(rep("black",2)))
abline(v=0,lty=3,col="gray50")

text(x=-2.2,y=1.5,labels="Linear fit",xpd=T,cex=1.3, adj = 0)
text(x=-2.2,y=0.5,labels="Polynomial fit",xpd=T,cex=1.3, adj = 0)
text(x=1,y= c(0.5,1.5),
     labels=c("[202; 1258]","[202; 1258]"),
     font=c(1),cex=1.2,xpd=T)

##################################################################
##################################################################
## Table S3.1: Numbers of papers and studies#####################

#avi.ES <- avi.ES[-which (avi.ES$Study_type == "Spatial"),] #excludes spatial studies to double check that results still hold throughout.  They do.

avi.animal <- avi.ES[which (avi.ES$Inv_kingdom == "Animal"),] 
avi.plant <- avi.ES[which (avi.ES$Inv_kingdom == "Plant"),] 

planthab.art <- avi.plant %>% group_by(Inv_habitat) %>% summarize(Art_count=length(unique(Article_ID)))
planthab.stud <- avi.plant %>% group_by(Inv_habitat) %>% summarize(Art_count=length(unique(Art_stud)))

anhab.art <- avi.animal %>% group_by(Inv_habitat) %>% summarize(Art_count=length(unique(Article_ID)))
anhab.stud <- avi.animal %>% group_by(Inv_habitat) %>% summarize(Art_count=length(unique(Art_stud)))

anhab.art
anhab.stud
planthab.art
planthab.stud

##################################################################
##################################################################
## Analysis 2: Single-species v. multi-species response ##########
# Analyze Multi_spp_resp which describes whether a study was population level or community level
# Drop response types that are "Other" in Response_type variable. 

#**NOTE: p values commented in code will change slightly when code is re-run, but should be similar to reported values**

avi.ES2 <- avi.ES[-which(avi.ES$Response_type == "Other"),]
n1 <- length(unique(avi.ES2$Article_ID[which(avi.ES2$Multi_spp_resp == "MULTIPLE")]))		# n1 = 156 articles (community level)
n2 <- length(unique(avi.ES2$Article_ID[which(avi.ES2$Multi_spp_resp == "SINGLE")]))
n3 <- length(unique(avi.ES2$Art_stud[which(avi.ES2$Multi_spp_resp == "MULTIPLE")]))
n4 <- length(unique(avi.ES2$Art_stud[which(avi.ES2$Multi_spp_resp == "SINGLE")]))

m2a <- MCMCglmm(ES.linear ~ Multi_spp_resp, random = ~ Article_ID, mev = avi.ES2$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES2)
m2a.int <- MCMCglmm(ES.linear ~ 1, random = ~ Article_ID, mev = avi.ES2$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES2)
#summary(m2a)
#significant community level linear effect ES = -0.7498
#population level effect = -0.7498 + 0.1325 = -0.6173
#population level effect does not differ significantly from community level effect p = 0.316

#estimate predicted linear effect sizes into object res2a

res2a <- cbind(m2a$Sol[, "(Intercept)"], m2a$Sol[, "(Intercept)"] + m2a$Sol[, "Multi_spp_respSINGLE"])

res2a.tmp <- HPDinterval(mcmc(res2a))
hpm <- c(mean(mcmc(res2a)[,1]),mean(mcmc(res2a)[,2]))
p.val <- c(pMCMC(mcmc(res2a)[,1]), pMCMC(mcmc(res2a)[,2]))
res2a <- data.frame(hpm,res2a.tmp, p.val)
colnames(res2a) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res2a) <- c("Community.l", "Population.l")

#random effects meta-analysis, polynomial term
m2b <- MCMCglmm(ES.poly ~ Multi_spp_resp, random = ~ Article_ID, mev = avi.ES2$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES2)
m2b.int <- MCMCglmm(ES.poly ~ 1, random = ~ Article_ID, mev = avi.ES2$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES2)

#summary(m2b)
#community level polynomial effect =  0.06651, does not differ sig from 0
#population level effect = 0.06651 + 0.18726 = 0.25377
#population level differs marginally non-significantly from community level effect p = 0.054

#estimate predicted poly effect sizes, into object res2b

res2b <- cbind(m2b$Sol[, "(Intercept)"], m2b$Sol[, "(Intercept)"] + m2b$Sol[, "Multi_spp_respSINGLE"])

res2b.tmp <- HPDinterval(mcmc(res2b))
hpm <- c(mean(mcmc(res2b)[,1]),mean(mcmc(res2b)[,2]))
p.val <- c(pMCMC(mcmc(res2b)[,1]), pMCMC(mcmc(res2b)[,2]))
res2b <- data.frame(hpm,res2b.tmp, p.val)
colnames(res2b) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res2b) <- c("Community.p", "Population.p")



##  Forest Plot Figure 2A Single-species Response ###############################

#Combine linear and polynomial estimates into a single table
res2a_2b <- rbind (res2b[2,],res2a[2,], res2b[1,], res2a[1,])

##New plot
par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(res2a_2b[2,1],res2a_2b[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.0,1.0),ylim=c(-1.0,1.0), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,1.5), col=0)
axis(side=4,labels=F,at=c(-1.5,1.5), col=0)
e.bars.points2(res2a_2b[2,1], z1 = res2a_2b[2,3], z2 = res2a_2b[2,2], x = c(res2a_2b[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(res2a_2b[1,1], z1 = res2a_2b[1,3], z2 = res2a_2b[1,2], x = c(res2a_2b[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(-1,-1,labels=paste0("[",n2,";",n4,"]"), pos=4)

#Model summaries, these contain an intercept term, so second coefficient is a difference, not a posterior mean
summary(m2a) #linear - no sig diff between single and multi-species response in linear term

##  Forest Plot Figure 2C Multi-species Response ###############################

res2c_2d <- res2a_2b[3:4,]

##New plot
par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(res2c_2d[2,1],res2c_2d[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.0,1.0),ylim=c(-1.0,1.0), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,1.5), col=0)
axis(side=4,labels=F,at=c(-1.5,1.5), col=0)
e.bars.points2(res2c_2d[2,1], z1 = res2c_2d[2,3], z2 = res2c_2d[2,2], x = c(res2c_2d[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(res2c_2d[1,1], z1 = res2c_2d[1,3], z2 = res2c_2d[1,2], x = c(res2c_2d[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(-1,-1,labels=paste0("[",n1,";",n3,"]"), pos=4)

#Model summaries, these contain an intercept term, so second coefficient is a difference, not a posterior mean
summary(m2b) #poly - single species response has a more positive polynomial than multi-species response (marginally non-significant p=0.07)


##################################################################
##################################################################
## Analysis 3: Trophic interactions for single vs. multi spp. response#####################
### Above = the invader is at a higher trophic level than the native
### Intra = the native and invader are on the same trophic level
### Below = the invader is at a lower trophic level than the native

#How much single-species data do we have?

# RERUN
avi.single <- avi.ES[which(avi.ES$Multi_spp_resp == "SINGLE"),]
single.art <- avi.single %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Article_ID)))
single.stud <- avi.single %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Art_stud)))

#Invader at higher trophic
n5 <- single.art[1,2]
n6 <- single.stud[1,2]

#Invader at same trophic
n7 <- single.art[5,2]
n8 <- single.stud[5,2]

#Invader at lower trophic
n15 <- single.art[3,2]
n16 <- single.stud[3,2]

#How much multi-species data do we have?
avi.multi <- avi.ES[which(avi.ES$Multi_spp_resp == "MULTIPLE"),]
multi.art <- avi.multi %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Article_ID)))
multi.stud <- avi.multi %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Art_stud)))

#Invader at higher trophic
n9 <- multi.art[1,2]
n10 <- multi.stud[1,2]

#Invader at same trophic
n11 <- multi.art[4,2]
n12 <- multi.stud[4,2]

#Invader at lower trophic
n13 <- multi.art[3,2]
n14 <- multi.stud[3,2]

# remove unwanted levels of Trophic_level
avi.ES3 <- avi.ES2[-c(which(avi.ES2$Trophic_level=="Above/Intra"),which(avi.ES2$Trophic_level=="Below/Intra"),which(avi.ES2$Trophic_level=="Mixed"),which(avi.ES2$Trophic_level=="Unknown")),]
# remove empty levels (clean up)
avi.ES3 <- droplevels (avi.ES3)

#model m7a
m7a <- MCMCglmm(ES.linear ~ Trophic_level + Multi_spp_resp, random = ~ Article_ID, mev = avi.ES3$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3)
summary(m7a)

# above and below sig diff in linear effect sizes p < 0.001
# above and intra sig diff in linear effect size p = 0.01

# set reference level of Trophic_level to intra, to examine intra vs. below
avi.ES3$Trophic_level <- relevel (avi.ES3$Trophic_level, ref = 3)

#model m7a2
m7a2 <- MCMCglmm(ES.linear ~ Trophic_level + Multi_spp_resp, random = ~ Article_ID, mev = avi.ES3$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3)
summary(m7a2)

# intra and below sig differ in linear effect sizes p < 0.001

# VERY IMPORTANT: Reset the Trophic_level variable to original indexing
avi.ES3$Trophic_level <- relevel (avi.ES3$Trophic_level, ref = 2)

#model m7b
m7b <- MCMCglmm(ES.poly ~ Trophic_level + Multi_spp_resp, random = ~ Article_ID, mev = avi.ES3$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3)
summary(m7b)

# above has significantly positive polynomial p = 0.002
# intra and below have significantly lower polynomial than above, p == 0.006 (intra), p = 0.004 (below)

# set reference level of Trophic_level to intra, to examine intra vs. below
avi.ES3$Trophic_level <- relevel (avi.ES3$Trophic_level, ref = 3)

#model m7b2
m7b2 <- MCMCglmm(ES.poly ~ Trophic_level + Multi_spp_resp, random = ~ Article_ID, mev = avi.ES3$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3)
summary(m7b2)

# intra and below do not sig differ in linear effect sizes p = 0.4

# VERY IMPORTANT: Reset the Trophic_level variable to original indexing
avi.ES3$Trophic_level <- relevel (avi.ES3$Trophic_level, ref = 2)

#estimate predicted linear effect sizes into object res7a

res7a <- cbind(m7a$Sol[, "(Intercept)"],
m7a$Sol[, "(Intercept)"] + m7a$Sol[, "Trophic_levelIntra"],
m7a$Sol[, "(Intercept)"] + m7a$Sol[, "Trophic_levelBelow"],
m7a$Sol[, "(Intercept)"] + m7a$Sol[, "Multi_spp_respSINGLE"],
m7a$Sol[, "(Intercept)"] + m7a$Sol[, "Trophic_levelIntra"]+ m7a$Sol[, "Multi_spp_respSINGLE"],
m7a$Sol[, "(Intercept)"] + m7a$Sol[, "Trophic_levelBelow"]+ m7a$Sol[, "Multi_spp_respSINGLE"])

res7a.tmp <- HPDinterval(mcmc(res7a))
hpm <- c(mean(mcmc(res7a)[,1]),mean(mcmc(res7a)[,2]),mean(mcmc(res7a)[,3]),mean(mcmc(res7a)[,4]),mean(mcmc(res7a)[,5]),mean(mcmc(res7a)[,6]))
p.val <- c(pMCMC(mcmc(res7a)[,1]), pMCMC(mcmc(res7a)[,2]), pMCMC(mcmc(res7a)[,3]), pMCMC(mcmc(res7a)[,4]), pMCMC(mcmc(res7a)[,5]), pMCMC(mcmc(res7a)[,6]))
res7a <- data.frame(hpm,res7a.tmp, p.val)
colnames(res7a) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res7a) <- c("Community.above.l","Community.intra.l","Community.below.l", "Popn.above.l","Popn.intra.l","Popn.below.l")

#estimate predicted polynomial effect sizes into object res7b

# First collect together all the posteriors we need, six in total, these are calculated from combinations of the posteriors corresponding to rows in the model summary
res7b <- cbind(m7b$Sol[, "(Intercept)"],
m7b$Sol[, "(Intercept)"] + m7b$Sol[, "Trophic_levelIntra"],
m7b$Sol[, "(Intercept)"] + m7b$Sol[, "Trophic_levelBelow"],
m7b$Sol[, "(Intercept)"] + m7b$Sol[, "Multi_spp_respSINGLE"],
m7b$Sol[, "(Intercept)"] + m7b$Sol[, "Trophic_levelIntra"]+ m7b$Sol[, "Multi_spp_respSINGLE"],
m7b$Sol[, "(Intercept)"] + m7b$Sol[, "Trophic_levelBelow"]+ m7b$Sol[, "Multi_spp_respSINGLE"])

res7b.tmp <- HPDinterval(mcmc(res7b))
hpm <- c(mean(mcmc(res7b)[,1]),mean(mcmc(res7b)[,2]),mean(mcmc(res7b)[,3]),mean(mcmc(res7b)[,4]),mean(mcmc(res7b)[,5]),mean(mcmc(res7b)[,6]))
p.val <- c(pMCMC(mcmc(res7b)[,1]), pMCMC(mcmc(res7b)[,2]), pMCMC(mcmc(res7b)[,3]), pMCMC(mcmc(res7b)[,4]), pMCMC(mcmc(res7b)[,5]), pMCMC(mcmc(res7b)[,6]))
res7b <- data.frame(hpm,res7b.tmp, p.val)
colnames(res7b) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res7b) <- c("Community.above.p","Community.intra.p","Community.below.p", "Popn.above.p","Popn.intra.p","Popn.below.p")



###############################

res7A <- rbind (res7b[4,],res7a[4,])

#Figure 3A invader at higher trophic
par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(res7A[2,1],res7A[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(res7A[2,1], z1 = res7A[2,3], z2 = res7A[2,2], x = c(res7A[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(res7A[1,1], z1 = res7A[1,3], z2 = res7A[1,2], x = c(res7A[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n5,";",n6,"]"), pos=4)

res7B <- rbind (res7b[5,],res7a[5,])

#Figure 3B invader at same trophic
par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(res7B[2,1],res7B[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(res7B[2,1], z1 = res7B[2,3], z2 = res7B[2,2], x = c(res7B[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(res7B[1,1], z1 = res7B[1,3], z2 = res7B[1,2], x = c(res7B[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n7,";",n8,"]"), pos=4)

res7C <- rbind (res7b[6,],res7a[6,])

#Figure 3C invader at lower trophic
par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(res7C[2,1],res7C[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(res7C[2,1], z1 = res7C[2,3], z2 = res7C[2,2], x = c(res7C[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(res7C[1,1], z1 = res7C[1,3], z2 = res7C[1,2], x = c(res7C[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n15,";",n16,"]"), pos=4)

res7D <- rbind (res7b[1,],res7a[1,])

#Figure 3D invader at higher trophic, multi species response
par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(res7D[2,1],res7D[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(res7D[2,1], z1 = res7D[2,3], z2 = res7D[2,2], x = c(res7D[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(res7D[1,1], z1 = res7D[1,3], z2 = res7D[1,2], x = c(res7D[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n9,";",n10,"]"), pos=4)

res7E <- rbind (res7b[2,],res7a[2,])

#Figure 3E invader at same trophic, multi species response
par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(res7E[2,1],res7E[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(res7E[2,1], z1 = res7E[2,3], z2 = res7E[2,2], x = c(res7E[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(res7E[1,1], z1 = res7E[1,3], z2 = res7E[1,2], x = c(res7E[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n11,";",n12,"]"), pos=4)

res7F <- rbind (res7b[3,],res7b[3,])

#Figure 3F invader at lower trophic, multi species response
par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(res7F[2,1],res7F[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(res7F[2,1], z1 = res7F[2,3], z2 = res7F[2,2], x = c(res7F[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(res7F[1,1], z1 = res7F[1,3], z2 = res7F[1,2], x = c(res7F[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n13,";",n14,"]"), pos=4)

##################################################################
##################################################################
## Analysis 4: Multi-species analysis by response type (Diversity, Evenness, Richness)
# not shown in manuscript

avi.ES4 <- avi.ES[which(avi.ES$Multi_spp_resp == "MULTIPLE"),]
avi.ES4 <- avi.ES4[-c(which(avi.ES4$Response_type=="Abundance"),which(avi.ES4$Response_type=="Other")),]
# remove empty levels (clean up)
avi.ES4 <- droplevels (avi.ES4)

# random effects meta-analysis, linear term, Response_type as fixed effects.
m4a <- MCMCglmm(ES.linear ~ Response_type, random = ~ Article_ID, mev = avi.ES4$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES4)
summary(m4a)

# Diversity has significantly negative linear effect size (P < 0.001)
# Richness and evenness do not differ significantly from diversity (linear effect) p = 0.308 (evenness), p = 0.066 (richness; marginally non-sig)

# set reference level of Response_type to richness, to examine richness vs. evenness
avi.ES4$Response_type <- relevel (avi.ES4$Response_type, ref = 3)

#model m4a2
m4a2 <- MCMCglmm(ES.linear ~ Response_type, random = ~ Article_ID, mev = avi.ES4$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES4)
summary(m4a2)

# richness and evenness sig differ in linear effect sizes p = 0.002

# VERY IMPORTANT: Reset the Response_type variable to original indexing
avi.ES4$Response_type <- relevel (avi.ES4$Response_type, ref = 2)


#random effects meta-analysis, polynomial term
m4b <- MCMCglmm(ES.poly ~ Response_type, random = ~ Article_ID, mev = avi.ES4$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES4)
summary(m4b)

# diversity polynomial effect does not differ significantly from zero (p = 0.514)
# evenness and richness do not differ from diversity in linear effect (p = 0.184, 0.564 respectively)

# set reference level of Response_type to richness, to examine richness vs. evenness
avi.ES4$Response_type <- relevel (avi.ES4$Response_type, ref = 3)

#model m4b2
m4b2 <- MCMCglmm(ES.poly ~ Response_type, random = ~ Article_ID, mev = avi.ES4$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES4)
summary(m4b2)

# Richness does not have a significant positive polynomial (p=0.11)
# richness and evenness  differ marginally non-significantly in polynomial effect sizes p = 0.06

# VERY IMPORTANT: Reset the Response_type variable to original indexing
avi.ES4$Response_type <- relevel (avi.ES4$Response_type, ref = 2)

#estimate predicted linear effect sizes into object res4a

res4a <- cbind(m4a$Sol[, "(Intercept)"] + m4a$Sol[, "Response_typeRichness"],
m4a$Sol[, "(Intercept)"], 
m4a$Sol[, "(Intercept)"]+ m4a$Sol[, "Response_typeEvenness"])

res4a.tmp <- HPDinterval(mcmc(res4a))
hpm <- c(mean(mcmc(res4a)[,1]),mean(mcmc(res4a)[,2]),mean(mcmc(res4a)[,3]))
p.val <- c(pMCMC(mcmc(res4a)[,1]), pMCMC(mcmc(res4a)[,2]), pMCMC(mcmc(res4a)[,3]))
res4a <- data.frame(hpm,res4a.tmp, p.val)
colnames(res4a) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res4a) <- c("Richness.l", "Diversity.l", "Evenness.l")


#estimate predicted polynomialr effect sizes into object res4b

res4b <- cbind(m4b$Sol[, "(Intercept)"] + m4b$Sol[, "Response_typeRichness"],
m4b$Sol[, "(Intercept)"], 
m4b$Sol[, "(Intercept)"] + m4b$Sol[, "Response_typeEvenness"])

res4b.tmp <- HPDinterval(mcmc(res4b))
hpm <- c(mean(mcmc(res4b)[,1]),mean(mcmc(res4b)[,2]),mean(mcmc(res4b)[,3]))
p.val <- c(pMCMC(mcmc(res4b)[,1]), pMCMC(mcmc(res4b)[,2]), pMCMC(mcmc(res4b)[,3]))
res4b <- data.frame(hpm,res4b.tmp, p.val)
colnames(res4b) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res4b) <- c("Richness.l", "Diversity.l", "Evenness.l")



##  Forest Plot Figure 4 Multi-spp by Response Type ###############################
#sample size:
commres.art <- avi.ES4 %>% group_by(Response_type) %>% summarize(Art_count=length(unique(Article_ID)))
commres.stud <- avi.ES4 %>% group_by(Response_type) %>% summarize(Art_count=length(unique(Art_stud)))

even.art <- commres.art[1,2]
even.stud <- commres.stud[1,2]

div.art <- commres.art[3,2]
div.stud <- commres.stud[3,2]

rich.art <- commres.art[2,2]
rich.stud <- commres.stud[2,2]

#Figure 4 linear, polynomial effects of richness, diversity & evenness
res4rich <- rbind (res4b[1,],res4a[1,])
res4div <- rbind (res4b[2,],res4a[2,])
res4even <- rbind (res4b[3,],res4a[3,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square
plot(res4rich[2,1],res4rich[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(res4rich[2,1], z1 = res4rich[2,3], z2 = res4rich[2,2], x = c(res4rich[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(res4rich[1,1], z1 = res4rich[1,3], z2 = res4rich[1,2], x = c(res4rich[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(-1,1,labels=paste0("[Richness=",rich.art,";",rich.stud,"]"), pos=4)
points(res4div[2,1],res4div[1,1],pch=17,col="red",cex=1.5)
e.bars.points2(res4div[2,1], z1 = res4div[2,3], z2 = res4div[2,2], x = c(res4div[1,1]),cap=0.1,clr=c(rep("red",2)), lwd=1.5)
e.bars.points3(res4div[1,1], z1 = res4div[1,3], z2 = res4div[1,2], x = c(res4div[2,1]),cap=0.1,clr=c(rep("red",2)), lwd=1.5, lower = T)
text(-1,.75,labels=paste0("[Diversity=",div.art,";",div.stud,"]"), pos=4)
points(res4even[2,1],res4even[1,1],pch=19,col="blue",cex=1.5)
e.bars.points2(res4even[2,1], z1 = res4even[2,3], z2 = res4even[2,2], x = c(res4even[1,1]),cap=0.1,clr=c(rep("blue",2)), lwd=1.5)
e.bars.points3(res4even[1,1], z1 = res4even[1,3], z2 = res4even[1,2], x = c(res4even[2,1]),cap=0.1,clr=c(rep("blue",2)), lwd=1.5, lower = T)
text(-1,.5,labels=paste0("[Evenness=",even.art,";",even.stud,"]"), pos=4)

summary(m4a)
summary(m4b)

##################################################################
##################################################################
## Analysis S3.3: Trophic interactions by recipient habitat#####################

avi.terr <- avi.ES[which (avi.ES$Inv_habitat == "terrestrial"),] 
avi.aqua <- avi.ES[which (avi.ES$Inv_habitat == "aquatic"),]
avi.marine <- avi.ES[which (avi.ES$Inv_habitat == "marine"),]

terr.art <- avi.terr %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Article_ID)))
terr.stud <- avi.terr %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Art_stud)))
#Invader at higher trophic
n17 <- terr.art[1,2]
n18 <- terr.stud[1,2]

#Invader at same trophic
n19 <- terr.art[5,2]
n20 <- terr.stud[5,2]

#Invader at lower trophic
n21 <- terr.art[3,2]
n22 <- terr.stud[3,2]

aqua.art <- avi.aqua %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Article_ID)))
aqua.stud <- avi.aqua %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Art_stud)))
#Invader at higher trophic
n23 <- aqua.art[1,2]
n24 <- aqua.stud[1,2]

#Invader at same trophic
n25 <- aqua.art[5,2]
n26 <- aqua.stud[5,2]

#Invader at lower trophic
n27 <- aqua.art[3,2]
n28 <- aqua.stud[3,2]

mar.art <- avi.marine %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Article_ID)))
mar.stud <- avi.marine %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Art_stud)))
#Invader at higher trophic
n29 <- mar.art[1,2]
n30 <- mar.stud[1,2]

#Invader at same trophic
n31 <- mar.art[4,2]
n32 <- mar.stud[4,2]

#Invader at lower trophic
n33 <- mar.art[2,2]
n34 <- mar.stud[2,2]


m5a <- MCMCglmm(ES.linear ~ Inv_habitat + Trophic_level, random = ~ Article_ID, mev = avi.ES3$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3)
summary(m5a)

# Linear effect for aquatic significantly differs from zero (negative effect) p < 0.001
# Linear effects for terr and marine do not differ significantly from linear effect for aquatic (p = 0.152; p = 0.812)

# set reference level of Inv_habitat to terr, to examine terr vs. marine
avi.ES3$Inv_habitat <- relevel (avi.ES3$Inv_habitat, ref = 3)

#model m5a2
m5a2 <- MCMCglmm(ES.linear ~ Inv_habitat + Trophic_level, random = ~ Article_ID, mev = avi.ES3$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3)
summary(m5a2)

# there is a marginally non-significant difference in linear effect between terr and marine habitats p = 0.086

# VERY IMPORTANT: Reset the Inv_habitat variable to original indexing
avi.ES3$Inv_habitat <- relevel (avi.ES3$Inv_habitat, ref = 2)



m5b <- MCMCglmm(ES.poly ~ Inv_habitat + Trophic_level, random = ~ Article_ID, mev = avi.ES3$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3)
summary(m5b)

# Poly effect for aquatic significantly differs from zero (positive effect) p < 0.001
# Poly effects for terr and marine do not differ significantly from poly effect for aquatic (p = 0.074; p = 0.250)

# set reference level of Inv_habitat to terr, to examine terr vs. marine
avi.ES3$Inv_habitat <- relevel (avi.ES3$Inv_habitat, ref = 3)

#model m5b2
m5b2 <- MCMCglmm(ES.poly ~ Inv_habitat + Trophic_level, random = ~ Article_ID, mev = avi.ES3$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES3)
summary(m5b2)

# there is no significant difference in polynomial effect between terr and marine habitats p = 0.826

# VERY IMPORTANT: Reset the Inv_habitat variable to original indexing
avi.ES3$Inv_habitat <- relevel (avi.ES3$Inv_habitat, ref = 2)



#estimate predicted linear effect sizes, into object res5a

# First collect together all the posteriors we need, nine in total, these are calculated from combinations of the posteriors corresponding to rows in the model summary
res5a <- cbind(m5a$Sol[, "(Intercept)"] + m5a$Sol[, "Inv_habitatterrestrial"],
m5a$Sol[, "(Intercept)"] + m5a$Sol[, "Inv_habitatterrestrial"]+ m5a$Sol[, "Trophic_levelIntra"],
m5a$Sol[, "(Intercept)"] + m5a$Sol[, "Inv_habitatterrestrial"]+ m5a$Sol[, "Trophic_levelBelow"],
m5a$Sol[, "(Intercept)"],
m5a$Sol[, "(Intercept)"] + m5a$Sol[, "Trophic_levelIntra"],
m5a$Sol[, "(Intercept)"] + m5a$Sol[, "Trophic_levelBelow"],
m5a$Sol[, "(Intercept)"]+ m5a$Sol[, "Inv_habitatmarine"],
m5a$Sol[, "(Intercept)"]+ m5a$Sol[, "Inv_habitatmarine"]+ m5a$Sol[, "Trophic_levelIntra"],
m5a$Sol[, "(Intercept)"]+ m5a$Sol[, "Inv_habitatmarine"]+ m5a$Sol[, "Trophic_levelBelow"])


res5a.tmp <- HPDinterval(mcmc(res5a))
hpm <- c(mean(mcmc(res5a)[,1]),mean(mcmc(res5a)[,2]),mean(mcmc(res5a)[,3]),mean(mcmc(res5a)[,4]),mean(mcmc(res5a)[,5]),mean(mcmc(res5a)[,6]),mean(mcmc(res5a)[,7]),mean(mcmc(res5a)[,8]),mean(mcmc(res5a)[,9]))
p.val <- c(pMCMC(mcmc(res5a)[,1]), pMCMC(mcmc(res5a)[,2]), pMCMC(mcmc(res5a)[,3]), pMCMC(mcmc(res5a)[,4]), pMCMC(mcmc(res5a)[,5]), pMCMC(mcmc(res5a)[,6]), pMCMC(mcmc(res5a)[,7]), pMCMC(mcmc(res5a)[,8]), pMCMC(mcmc(res5a)[,9]))
res5a <- data.frame(hpm,res5a.tmp, p.val)
colnames(res5a) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res5a) <- c("Terr.above.l","Terr.intra.l","Terr.below.l","Fresh.above.l","Fresh.intra.l","Fresh.below.l", "Mar.above.l","Mar.intra.l","Mar.below.l")

#estimate predicted polynomial effect sizes, into object res5b

# First collect together all the posteriors we need, nine in total, these are calculated from combinations of the posteriors corresponding to rows in the model summary
res5b <- cbind(m5b$Sol[, "(Intercept)"] + m5b$Sol[, "Inv_habitatterrestrial"],
m5b$Sol[, "(Intercept)"] + m5b$Sol[, "Inv_habitatterrestrial"]+ m5b$Sol[, "Trophic_levelIntra"],
m5b$Sol[, "(Intercept)"] + m5b$Sol[, "Inv_habitatterrestrial"]+ m5b$Sol[, "Trophic_levelBelow"],
m5b$Sol[, "(Intercept)"],
m5b$Sol[, "(Intercept)"] + m5b$Sol[, "Trophic_levelIntra"],
m5b$Sol[, "(Intercept)"] + m5b$Sol[, "Trophic_levelBelow"],
m5b$Sol[, "(Intercept)"]+ m5b$Sol[, "Inv_habitatmarine"],
m5b$Sol[, "(Intercept)"]+ m5b$Sol[, "Inv_habitatmarine"]+ m5b$Sol[, "Trophic_levelIntra"],
m5b$Sol[, "(Intercept)"]+ m5b$Sol[, "Inv_habitatmarine"]+ m5b$Sol[, "Trophic_levelBelow"])


res5b.tmp <- HPDinterval(mcmc(res5b))
hpm <- c(mean(mcmc(res5b)[,1]),mean(mcmc(res5b)[,2]),mean(mcmc(res5b)[,3]),mean(mcmc(res5b)[,4]),mean(mcmc(res5b)[,5]),mean(mcmc(res5b)[,6]),mean(mcmc(res5b)[,7]),mean(mcmc(res5b)[,8]),mean(mcmc(res5b)[,9]))
p.val <- c(pMCMC(mcmc(res5b)[,1]), pMCMC(mcmc(res5b)[,2]), pMCMC(mcmc(res5b)[,3]), pMCMC(mcmc(res5b)[,4]), pMCMC(mcmc(res5b)[,5]), pMCMC(mcmc(res5b)[,6]), pMCMC(mcmc(res5b)[,7]), pMCMC(mcmc(res5b)[,8]), pMCMC(mcmc(res5b)[,9]))
res5b <- data.frame(hpm,res5b.tmp, p.val)
colnames(res5b) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res5b) <- c("Terr.above.p","Terr.intra.p","Terr.below.p","Fresh.above.p","Fresh.intra.p","Fresh.below.p", "Mar.above.p","Mar.intra.p","Mar.below.p")




## Forest Plot Figure S3.3a-f Trophic interactions by recipient habitat ###############################

#Figure S3.3A terrestrial invader at higher trophic
resS2.2A <- rbind (res5b[1,],res5a[1,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.2A[2,1],resS2.2A[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.2A[2,1], z1 = resS2.2A[2,3], z2 = resS2.2A[2,2], x = c(resS2.2A[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.2A[1,1], z1 = resS2.2A[1,3], z2 = resS2.2A[1,2], x = c(resS2.2A[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n17,";",n18,"]"), pos=4)

#Figure S3.3B terrestrial invader at same trophic
resS2.2B <- rbind (res5b[2,],res5a[2,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.2B[2,1],resS2.2B[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.2B[2,1], z1 = resS2.2B[2,3], z2 = resS2.2B[2,2], x = c(resS2.2B[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.2B[1,1], z1 = resS2.2B[1,3], z2 = resS2.2B[1,2], x = c(resS2.2B[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n19,";",n20,"]"), pos=4)

#Figure S3.3C terrestrial invader at lower trophic
resS2.2C <- rbind (res5b[3,],res5a[3,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.2C[2,1],resS2.2C[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.2C[2,1], z1 = resS2.2C[2,3], z2 = resS2.2C[2,2], x = c(resS2.2C[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.2C[1,1], z1 = resS2.2C[1,3], z2 = resS2.2C[1,2], x = c(resS2.2C[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n21,";",n22,"]"), pos=4)

#Figure S3.3D aquatic invader at higher trophic
resS2.2D <- rbind (res5b[4,],res5a[4,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.2D[2,1],resS2.2D[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.2D[2,1], z1 = resS2.2D[2,3], z2 = resS2.2D[2,2], x = c(resS2.2D[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.2D[1,1], z1 = resS2.2D[1,3], z2 = resS2.2D[1,2], x = c(resS2.2D[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n23,";",n24,"]"), pos=4)

#Figure S3.3E aquatic invader at same trophic
resS2.2E <- rbind (res5b[5,],res5a[5,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.2E[2,1],resS2.2E[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.2E[2,1], z1 = resS2.2E[2,3], z2 = resS2.2E[2,2], x = c(resS2.2E[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.2E[1,1], z1 = resS2.2E[1,3], z2 = resS2.2E[1,2], x = c(resS2.2E[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n25,";",n26,"]"), pos=4)

#Figure S3.3F aquatic invader at lower trophic
resS2.2F <- rbind (res5b[6,],res5a[6,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.2F[2,1],resS2.2F[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.2F[2,1], z1 = resS2.2F[2,3], z2 = resS2.2F[2,2], x = c(resS2.2F[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.2F[1,1], z1 = resS2.2F[1,3], z2 = resS2.2F[1,2], x = c(resS2.2F[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n27,";",n28,"]"), pos=4)

#Figure S3.3G marine invader at higher trophic
resS2.2G <- rbind (res5b[7,],res5a[7,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.2G[2,1],resS2.2G[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.2G[2,1], z1 = resS2.2G[2,3], z2 = resS2.2G[2,2], x = c(resS2.2G[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.2G[1,1], z1 = resS2.2G[1,3], z2 = resS2.2G[1,2], x = c(resS2.2G[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n29,";",n30,"]"), pos=4)

#Figure S3.3H marine invader at same trophic
resS2.2H <- rbind (res5b[8,],res5a[8,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.2H[2,1],resS2.2H[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.2H[2,1], z1 = resS2.2H[2,3], z2 = resS2.2H[2,2], x = c(resS2.2H[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.2H[1,1], z1 = resS2.2H[1,3], z2 = resS2.2H[1,2], x = c(resS2.2H[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n31,";",n32,"]"), pos=4)

#Figure S3.3I marine invader at lower trophic
resS2.2I <- rbind (res5b[9,],res5a[9,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.2I[2,1],resS2.2I[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.2I[2,1], z1 = resS2.2I[2,3], z2 = resS2.2I[2,2], x = c(resS2.2I[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.2I[1,1], z1 = resS2.2I[1,3], z2 = resS2.2I[1,2], x = c(resS2.2I[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n33,";",n34,"]"), pos=4)

res5a
res5b

##################################################################
##################################################################
## Analysis S3.4: Trophic interactions by invader taxon#####################


avi.ES5 <- avi.ES3[-which(avi.ES3$Inv_kingdom =="Bacteria"),]
avi.ES5 <- droplevels(avi.ES5)

avi.animal <- avi.ES5[which(avi.ES5$Inv_kingdom =="Animal"),]
an.art <- avi.animal %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Article_ID)))
an.stud <- avi.animal %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Art_stud)))
#Invader at higher trophic
n35 <- an.art[1,2]
n36 <- an.stud[1,2]

#Invader at same trophic
n37 <- an.art[3,2]
n38 <- an.stud[3,2]

#Invader at lower trophic
n39 <- an.art[2,2]
n40 <- an.stud[2,2]

avi.plant <- avi.ES5[which(avi.ES5$Inv_kingdom =="Plant"),]
plant.art <- avi.plant %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Article_ID)))
plant.stud <- avi.plant %>% group_by(Trophic_level) %>% summarize(Art_count=length(unique(Art_stud)))
#Invader at same trophic
n41 <- plant.art[2,2]
n42 <- plant.stud[2,2]

#Invader at lower trophic
n43 <- plant.art[1,2]
n44 <- plant.stud[1,2]

m6a <- MCMCglmm(ES.linear ~ Inv_kingdom + Trophic_level, random = ~ Article_ID, mev = avi.ES5$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES5)
summary(m6a)

# Linear effect size for invasive animals is significantly less than zero (p < 0.001)
# Linear effect sizes significantly less negative in invasive plants than invasive animals (p = 0.002)

m6b <- MCMCglmm(ES.poly ~ Inv_kingdom + Trophic_level, random = ~ Article_ID, mev = avi.ES5$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES5)
summary(m6b)

# Polynomial effect size for invasive animals is significantly greater than zero (p < 0.001)
# Polynomial effect sizes do not differ between invasive plants and invasive animals (p = 0.104)


#estimate predicted linear effect sizes, into object res6a

res6a <- cbind(m6a$Sol[, "(Intercept)"],
m6a$Sol[, "(Intercept)"] + m6a$Sol[, "Trophic_levelIntra"],
m6a$Sol[, "(Intercept)"] + m6a$Sol[, "Trophic_levelBelow"],
m6a$Sol[, "(Intercept)"] + m6a$Sol[, "Inv_kingdomPlant"],
m6a$Sol[, "(Intercept)"] + m6a$Sol[, "Trophic_levelIntra"]+ m6a$Sol[, "Inv_kingdomPlant"],
m6a$Sol[, "(Intercept)"] + m6a$Sol[, "Trophic_levelBelow"]+ m6a$Sol[, "Inv_kingdomPlant"])

res6a.tmp <- HPDinterval(mcmc(res6a))
hpm <- c(mean(mcmc(res6a)[,1]),mean(mcmc(res6a)[,2]),mean(mcmc(res6a)[,3]),mean(mcmc(res6a)[,4]),mean(mcmc(res6a)[,5]),mean(mcmc(res6a)[,6]))
p.val <- c(pMCMC(mcmc(res6a)[,1]), pMCMC(mcmc(res6a)[,2]), pMCMC(mcmc(res6a)[,3]), pMCMC(mcmc(res6a)[,4]), pMCMC(mcmc(res6a)[,5]), pMCMC(mcmc(res6a)[,6]))
res6a <- data.frame(hpm,res6a.tmp, p.val)
colnames(res6a) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res6a) <- c("Animal.above.l","Animal.intra.l","Animal.below.l", "Plant.above.l","Plant.intra.l","Plant.below.l")

#estimate predicted polynomial effect sizes into object res6b

# First collect together all the posteriors we need, six in total, these are calculated from combinations of the posteriors corresponding to rows in the model summary
res6b <- cbind(m6b$Sol[, "(Intercept)"],
m6b$Sol[, "(Intercept)"] + m6b$Sol[, "Trophic_levelIntra"],
m6b$Sol[, "(Intercept)"] + m6b$Sol[, "Trophic_levelBelow"],
m6b$Sol[, "(Intercept)"] + m6b$Sol[, "Inv_kingdomPlant"],
m6b$Sol[, "(Intercept)"] + m6b$Sol[, "Trophic_levelIntra"]+ m6b$Sol[, "Inv_kingdomPlant"],
m6b$Sol[, "(Intercept)"] + m6b$Sol[, "Trophic_levelBelow"]+ m6b$Sol[, "Inv_kingdomPlant"])

res6b.tmp <- HPDinterval(mcmc(res6b))
hpm <- c(mean(mcmc(res6b)[,1]),mean(mcmc(res6b)[,2]),mean(mcmc(res6b)[,3]),mean(mcmc(res6b)[,4]),mean(mcmc(res6b)[,5]),mean(mcmc(res6b)[,6]))
p.val <- c(pMCMC(mcmc(res6b)[,1]), pMCMC(mcmc(res6b)[,2]), pMCMC(mcmc(res6b)[,3]), pMCMC(mcmc(res6b)[,4]), pMCMC(mcmc(res6b)[,5]), pMCMC(mcmc(res6b)[,6]))
res6b <- data.frame(hpm,res6b.tmp, p.val)
colnames(res6b) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res6b) <- c("Animal.above.p","Animal.intra.p","Animal.below.p", "Plant.above.p","Plant.intra.p","Plant.below.p")


#### Forest Plot Figure S3.4a-d Trophic interactions by invader taxon ###############################

#Figure S3.4A animal invader at higher trophic
resS2.3A <- rbind (res6b[1,],res6a[1,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.3A[2,1],resS2.3A[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.3A[2,1], z1 = resS2.3A[2,3], z2 = resS2.3A[2,2], x = c(resS2.3A[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.3A[1,1], z1 = resS2.3A[1,3], z2 = resS2.3A[1,2], x = c(resS2.3A[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n35,";",n36,"]"), pos=4)

#Figure S3.4B animal invader at same trophic
resS2.3B <- rbind (res6b[2,],res6a[2,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.3B[2,1],resS2.3B[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.3B[2,1], z1 = resS2.3B[2,3], z2 = resS2.3B[2,2], x = c(resS2.3B[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.3B[1,1], z1 = resS2.3B[1,3], z2 = resS2.3B[1,2], x = c(resS2.3B[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n37,";",n38,"]"), pos=4)

#Figure S3.4C animal invader at lower trophic
resS2.3C <- rbind (res6b[3,],res6a[3,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.3C[2,1],resS2.3C[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.3C[2,1], z1 = resS2.3C[2,3], z2 = resS2.3C[2,2], x = c(resS2.3C[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.3C[1,1], z1 = resS2.3C[1,3], z2 = resS2.3C[1,2], x = c(resS2.3C[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n39,";",n40,"]"), pos=4)

#Figure S3.4D plant invader at same trophic
resS2.3D <- rbind (res6b[5,],res6a[5,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.3D[2,1],resS2.3D[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.3D[2,1], z1 = resS2.3D[2,3], z2 = resS2.3D[2,2], x = c(resS2.3D[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.3D[1,1], z1 = resS2.3D[1,3], z2 = resS2.3D[1,2], x = c(resS2.3D[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n41,";",n42,"]"), pos=4)

#Figure S3.4E plant invader at lower trophic
resS2.3E <- rbind (res6b[6,],res6a[6,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.3E[2,1],resS2.3E[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=22,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.3E[2,1], z1 = resS2.3E[2,3], z2 = resS2.3E[2,2], x = c(resS2.3E[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.3E[1,1], z1 = resS2.3E[1,3], z2 = resS2.3E[1,2], x = c(resS2.3E[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(0,1,labels=paste0("[",n43,";",n44,"]"), pos=4)

res6a #linear
res6b #poly

##################################################################
##################################################################
## Analysis S3.5A: Effect by study type#####################


avi.ES6 <- avi.ES[-which(avi.ES$Study_type=="Multiple"),]
avi.ES6 <- droplevels(avi.ES6)

m8a <- MCMCglmm(ES.linear ~ Study_type, random = ~ Article_ID, mev = avi.ES6$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES6)
summary(m8a)

#Experimental studies had significantly negative linear effect size (p < 0.001)
#Neither of the Linear effects for Spatial or temporal studies differed significantly from that for experimental studies (p = 0.696; p = 0.434)


# set reference level of Study_type to Spatial, to examine Spatial vs. Temporal
avi.ES6$Study_type <- relevel (avi.ES6$Study_type, ref = 2)

#model m8a2
m8a2 <- MCMCglmm(ES.linear ~ Study_type, random = ~ Article_ID, mev = avi.ES6$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES6)
summary(m8a2)

# there was no significant difference in polynomial effect between spatial and temporal studies p = 0.5

# VERY IMPORTANT: Reset the Study_type variable to original indexing
avi.ES6$Study_type <- relevel (avi.ES6$Study_type, ref = 2)



#random effects meta-analysis, polynomial term
m8b <- MCMCglmm(ES.poly ~ Study_type, random = ~ Article_ID, mev = avi.ES6$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES6)
summary(m8b)

# The polynomial effect for experimental studies was not significantly different from zero (p = 0.440)
# #Neither of the polynomial effects for Spatial or temporal studies differed significantly from that for experimental studies (p = 0.79; p = 0.11)

# set reference level of Study_type to Spatial, to examine Spatial vs. Temporal
avi.ES6$Study_type <- relevel (avi.ES6$Study_type, ref = 2)

#model m8b2
m8b2 <- MCMCglmm(ES.poly ~ Study_type, random = ~ Article_ID, mev = avi.ES6$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES6)
summary(m8b2)

# The polynomial effect differed significantly between spatial and temporal studies p = 0.02

# VERY IMPORTANT: Reset the Study_type variable to original indexing
avi.ES6$Study_type <- relevel (avi.ES6$Study_type, ref = 2)



#estimate predicted linear effect sizes into object res8a

res8a <- cbind(m8a$Sol[, "(Intercept)"], m8a$Sol[, "(Intercept)"] + m8a$Sol[, "Study_typeSpatial"], m8a$Sol[, "(Intercept)"] + m8a$Sol[, "Study_typeTemporal"])

res8a.tmp <- HPDinterval(mcmc(res8a))
hpm <- c(mean(mcmc(res8a)[,1]),mean(mcmc(res8a)[,2]),mean(mcmc(res8a)[,3]))
p.val <- c(pMCMC(mcmc(res8a)[,1]), pMCMC(mcmc(res8a)[,2]), pMCMC(mcmc(res8a)[,3]))
res8a <- data.frame(hpm,res8a.tmp, p.val)
colnames(res8a) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res8a) <- c("Exptal.l", "Spatial.l", "Temporal.p")



#estimate predicted poly effect sizes, into object res8b

res8b <- cbind(m8b$Sol[, "(Intercept)"], m8b$Sol[, "(Intercept)"] + m8b$Sol[, "Study_typeSpatial"], m8b$Sol[, "(Intercept)"] + m8b$Sol[, "Study_typeTemporal"])

res8b.tmp <- HPDinterval(mcmc(res8b))
hpm <- c(mean(mcmc(res8b)[,1]),mean(mcmc(res8b)[,2]),mean(mcmc(res8b)[,3]))
p.val <- c(pMCMC(mcmc(res8b)[,1]), pMCMC(mcmc(res8b)[,2]), pMCMC(mcmc(res8b)[,3]))
res8b <- data.frame(hpm,res8b.tmp, p.val)
colnames(res8b) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res8b) <- c("Exptal.p", "Spatial.p", "Temporal.p")



##  Bivariate Forest Plot Figure S3.5A Effect Size by Study Type ###############################
#sample size:
studytype.art <- avi.ES6 %>% group_by(Study_type) %>% summarize(Art_count=length(unique(Article_ID)))
studytype.stud <- avi.ES6 %>% group_by(Study_type) %>% summarize(Art_count=length(unique(Art_stud)))

exp.art <- studytype.art[1,2]
exp.stud <- studytype.stud[1,2]

spat.art <- studytype.art[2,2]
spat.stud <- studytype.stud[2,2]

temp.art <- studytype.art[3,2]
temp.stud <- studytype.stud[3,2]

resS2.4exp <- rbind (res8b[1,],res8a[1,])
resS2.4spat <- rbind (res8b[2,],res8a[2,])
resS2.4temp <- rbind (res8b[3,],res8a[3,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(resS2.4exp[2,1],resS2.4exp[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=17,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(resS2.4exp[2,1], z1 = resS2.4exp[2,3], z2 = resS2.4exp[2,2], x = c(resS2.4exp[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(resS2.4exp[1,1], z1 = resS2.4exp[1,3], z2 = resS2.4exp[1,2], x = c(resS2.4exp[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(-1.5,-1,labels=paste0("[Experimental=",exp.art,";",exp.stud,"]"), pos=4)
points(resS2.4spat[2,1],resS2.4spat[1,1],pch=15,col="red",cex=1.5)
e.bars.points2(resS2.4spat[2,1], z1 = resS2.4spat[2,3], z2 = resS2.4spat[2,2], x = c(resS2.4spat[1,1]),cap=0.1,clr=c(rep("red",2)), lwd=1.5)
e.bars.points3(resS2.4spat[1,1], z1 = resS2.4spat[1,3], z2 = resS2.4spat[1,2], x = c(resS2.4spat[2,1]),cap=0.1,clr=c(rep("red",2)), lwd=1.5, lower = T)
text(-1.5,-.75,labels=paste0("[Spatial=",spat.art,";",spat.stud,"]"), pos=4, col="red")
points(resS2.4temp[2,1],resS2.4temp[1,1],pch=19,col="blue",cex=1.5)
e.bars.points2(resS2.4temp[2,1], z1 = resS2.4temp[2,3], z2 = resS2.4temp[2,2], x = c(resS2.4temp[1,1]),cap=0.1,clr=c(rep("blue",2)), lwd=1.5)
e.bars.points3(resS2.4temp[1,1], z1 = resS2.4temp[1,3], z2 = resS2.4temp[1,2], x = c(resS2.4temp[2,1]),cap=0.1,clr=c(rep("blue",2)), lwd=1.5, lower = T)
text(-1.5,-.5,labels=paste0("[Temporal=",temp.art,";",temp.stud,"]"), pos=4, col="blue")

res8a #linear
res8b #poly

## Analysis S3.5B: Effect by study type#####################
#Analysis for only studies with known trophic interaction

avi.ES6B <- avi.ES3[-which(avi.ES3$Study_type=="Multiple"),]  #avi.ES3 are the data used for the trophic analyses
avi.ES6B <- droplevels(avi.ES6B)

m8ab <- MCMCglmm(ES.linear ~ Study_type, random = ~ Article_ID, mev = avi.ES6B$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES6B)
summary(m8ab)

#Experimental studies had significantly negative linear effect size (p < 0.001)
#Neither of the Linear effects for Spatial or temporal studies differed significantly from that for experimental studies (p = 0.58; p = 0.72)


# set reference level of Study_type to Spatial, to examine Spatial vs. Temporal
avi.ES6B$Study_type <- relevel (avi.ES6B$Study_type, ref = 2)

#model m8ab2
m8ab2 <- MCMCglmm(ES.linear ~ Study_type, random = ~ Article_ID, mev = avi.ES6B$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES6B)
summary(m8ab2)

# there was no significant difference in polynomial effect between spatial and temporal studies p = 0.35

# VERY IMPORTANT: Reset the Study_type variable to original indexing
avi.ES6B$Study_type <- relevel (avi.ES6B$Study_type, ref = 2)


#random effects meta-analysis, polynomial term
m8bb <- MCMCglmm(ES.poly ~ Study_type, random = ~ Article_ID, mev = avi.ES6B$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES6B)
summary(m8bb)

# The polynomial effect for experimental studies was not significantly different from zero (p = 0.480)
# #Neither of the polynomial effects for Spatial or temporal studies differed significantly from that for experimental studies (p = 0.93; p = 0.18)

# set reference level of Study_type to Spatial, to examine Spatial vs. Temporal
avi.ES6B$Study_type <- relevel (avi.ES6B$Study_type, ref = 2)

#model m8bb2
m8bb2 <- MCMCglmm(ES.poly ~ Study_type, random = ~ Article_ID, mev = avi.ES6B$mev, prior = prior, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES6B)
summary(m8bb2)

# The polynomial effect was not significantly different between spatial and temporal studies p = 0.17

# VERY IMPORTANT: Reset the Study_type variable to original indexing
avi.ES6B$Study_type <- relevel (avi.ES6B$Study_type, ref = 2)



#estimate predicted linear effect sizes into object res8ab

res8ab <- cbind(m8ab$Sol[, "(Intercept)"], m8ab$Sol[, "(Intercept)"] + m8ab$Sol[, "Study_typeSpatial"], m8ab$Sol[, "(Intercept)"] + m8ab$Sol[, "Study_typeTemporal"])

res8ab.tmp <- HPDinterval(mcmc(res8ab))
hpm <- c(mean(mcmc(res8ab)[,1]),mean(mcmc(res8ab)[,2]),mean(mcmc(res8ab)[,3]))
p.val <- c(pMCMC(mcmc(res8ab)[,1]), pMCMC(mcmc(res8ab)[,2]), pMCMC(mcmc(res8ab)[,3]))
res8ab <- data.frame(hpm,res8ab.tmp, p.val)
colnames(res8ab) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res8ab) <- c("Exptal.l", "Spatial.l", "Temporal.l")



#estimate predicted poly effect sizes, into object res8bb

res8bb <- cbind(m8bb$Sol[, "(Intercept)"], m8bb$Sol[, "(Intercept)"] + m8bb$Sol[, "Study_typeSpatial"], m8bb$Sol[, "(Intercept)"] + m8bb$Sol[, "Study_typeTemporal"])

res8bb.tmp <- HPDinterval(mcmc(res8bb))
hpm <- c(mean(mcmc(res8bb)[,1]),mean(mcmc(res8bb)[,2]),mean(mcmc(res8bb)[,3]))
p.val <- c(pMCMC(mcmc(res8bb)[,1]), pMCMC(mcmc(res8bb)[,2]), pMCMC(mcmc(res8bb)[,3]))
res8bb <- data.frame(hpm,res8bb.tmp, p.val)
colnames(res8bb) <- c("post.mean", "l-95% CI", "u-95% CI", "pMCMC")
rownames(res8bb) <- c("Exptal.p", "Spatial.p", "Temporal.p")



##  Bivariate Forest Plot Figure S3.5B Effect Size by Study Type ###############################
#sample size:
studytype.artb <- avi.ES6B %>% group_by(Study_type) %>% summarize(Art_count=length(unique(Article_ID)))
studytype.studb <- avi.ES6B %>% group_by(Study_type) %>% summarize(Art_count=length(unique(Art_stud)))

exp.artb <- studytype.artb[1,2]
exp.studb <- studytype.studb[1,2]

spat.artb <- studytype.artb[2,2]
spat.studb <- studytype.studb[2,2]

temp.artb <- studytype.artb[3,2]
temp.studb <- studytype.studb[3,2]

res2.4bexp <- rbind (res8bb[1,],res8ab[1,])
res2.4bspat <- rbind (res8bb[2,],res8ab[2,])
res2.4btemp <- rbind (res8bb[3,],res8ab[3,])

par(mfrow = c(1,1))
par(mar=c(5,4,2,2))
par(pty="s") #square!!
plot(res2.4bexp[2,1],res2.4bexp[1,1],xlab="Linear effect size", ylab="Polynomial effect size", bg=c(rep("black",2)), pch=17,
     cex=1.5, xlim=c(-1.5,0.5),ylim=c(-1,1), abline(h=0,lty=5,col="gray50", lwd=1.5), asp=1)
abline(v=0,lty=5,col="gray50", lwd=1.5)
axis(side=3,labels=F,at=c(-1.5,2), col=0)
axis(side=4,labels=F,at=c(-1.5,2), col=0)
e.bars.points2(res2.4bexp[2,1], z1 = res2.4bexp[2,3], z2 = res2.4bexp[2,2], x = c(res2.4bexp[1,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5)
e.bars.points3(res2.4bexp[1,1], z1 = res2.4bexp[1,3], z2 = res2.4bexp[1,2], x = c(res2.4bexp[2,1]),cap=0.1,clr=c(rep("black",2)), lwd=1.5, lower = T)
text(-1.5,-1,labels=paste0("[Experimental=",exp.artb,";",exp.studb,"]"), pos=4)
points(res2.4bspat[2,1],res2.4bspat[1,1],pch=15,col="red",cex=1.5)
e.bars.points2(res2.4bspat[2,1], z1 = res2.4bspat[2,3], z2 = res2.4bspat[2,2], x = c(res2.4bspat[1,1]),cap=0.1,clr=c(rep("red",2)), lwd=1.5)
e.bars.points3(res2.4bspat[1,1], z1 = res2.4bspat[1,3], z2 = res2.4bspat[1,2], x = c(res2.4bspat[2,1]),cap=0.1,clr=c(rep("red",2)), lwd=1.5, lower = T)
text(-1.5,-.75,labels=paste0("[Spatial=",spat.artb,";",spat.studb,"]"), pos=4, col="red")
points(res2.4btemp[2,1],res2.4btemp[1,1],pch=19,col="blue",cex=1.5)
e.bars.points2(res2.4btemp[2,1], z1 = res2.4btemp[2,3], z2 = res2.4btemp[2,2], x = c(res2.4btemp[1,1]),cap=0.1,clr=c(rep("blue",2)), lwd=1.5)
e.bars.points3(res2.4btemp[1,1], z1 = res2.4btemp[1,3], z2 = res2.4btemp[1,2], x = c(res2.4btemp[2,1]),cap=0.1,clr=c(rep("blue",2)), lwd=1.5, lower = T)
text(-1.5,-.5,labels=paste0("[Temporal=",temp.artb,";",temp.studb,"]"), pos=4, col="blue")

res8ab #linear
res8bb #poly

##################################################################
##################################################################
##################################################################
# SECTION 6. Run models for funnel plots

# Equivalent MCMCglmm code and metafor code for fitting random-effects meta-analysis on AvI data
#f1a <- MCMCglmm(ES.linear ~ 1, mev = avi.ES$mev, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES)
f1a <- rma (yi = ES.linear, vi = avi.ES$mev, mods = ~ 1, data = avi.ES, method = "REML")
summary(f1a)

#f1b <- MCMCglmm(ES.poly ~ 1, mev = avi.ES$mev, nitt = 110000, burnin = 10000, thin = 100, verbose = F, data = avi.ES)
f2a <- rma (yi = ES.poly, vi = avi.ES$mev, mods = ~ 1, data = avi.ES, method = "REML")
summary(f2a)

##################################################################
##################################################################
##################################################################
# SECTION 7. Funnel plots based on f1a, b above 

par (mfrow = c(1, 2))
funnel(f1a,level=c(90,95,99),shade=c("white", "gray", "darkgray"),addtau2=T, vi = avi.ES$mev,cex.axis=1.5,cex.lab=1.5,xlab="Linear effect size")
funnel(f2a,level=c(90,95,99),shade=c("white", "gray", "darkgray"),addtau2=T, vi = avi.ES$mev,cex.axis=1.5,cex.lab=1.5,xlab="Polynomial effect size")