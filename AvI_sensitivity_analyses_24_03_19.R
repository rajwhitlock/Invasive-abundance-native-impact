# This script presents additional sensitivity analyses that were devised and carried out by Raj Whitlock, that complement standard model checking for the Bayesian mixed effects analysis presented in files AvI_metaanalysis_partial-r_24_03_19.R and AvI_metaanalysis_slopes__24_03_19.R
##################################################################
##################################################################
##################################################################
# code to look see how the sign of the raw data is distributed, both for invasive and native species, run the hashed code following the function
loc.xy <- function (data) {
	
	xx <- split(data,list(data$Article_ID,data$Study_ID),drop=T)	
yy <- lapply(xx, function(y){
	c(sign(min(y$Response)), sign(max(y$Response)), sign(min(y$Abundance_Invader)), sign(max(y$Abundance_Invader)))
})

yy2 <- do.call(rbind,yy)

return (yy2)
	
}

# sign.data <- loc.xy(data)
# table(sign.data[,1],sign.data[,2])
# table(sign.data[,3],sign.data[,4])



##################################################################
##################################################################
##################################################################
# code to quantify the extent to which invasive species abundance data are distributed evenly over the range in invasive species abundance presented within each study

xx <- split(data,list(data$Article_ID,data$Study_ID),drop=T)	
yy <- lapply(xx, function(y){

if(sign (min (y$Abundance_Invader))==-1){offset.x <- min (y$Abundance_Invader)} else {offset.x <- 0}
		inv.recentred <- (y$Abundance_Invader - offset.x)/(max(y$Abundance_Invader) - offset.x)
		min(inv.recentred)/max(inv.recentred)})

yy <- do.call(rbind,yy)

hist(yy)

#What number of studies have a "full length" distribution of invasive abundance?
length(which(yy<0.2))
#1050
#1050/1115

plot(yy,samp.sizes, ylim=c(0,500))
# samp.sizes computed below

# code to quantify the extent to which invasive species abundance data are distributed evenly over the range in invasive species abundance presented within each study, using correspondence between median and mean as a metric

yy2 <- lapply(xx, function(y){

if(sign (min (y$Abundance_Invader))==-1){offset.x <- min (y$Abundance_Invader)} else {offset.x <- 0}
		inv.recentred <- (y$Abundance_Invader - offset.x)/(max(y$Abundance_Invader) - offset.x)
		log(median(inv.recentred)/mean(inv.recentred))})

yy2 <- do.call(rbind,yy2)

# how many studies have invasive abundance data with mean within 25% of the median value
length(which(-0.2876821 < yy2 & yy2 < 0.2231436))
#542

# how many studies have invasive abundance data with mean within 50% of the median value
length(which(-0.6931472 < yy2 & yy2 < 0.4054651))
#829
#829/1115

plot(yy2,samp.sizes, ylim=c(0,500))
# samp.sizes computed below

plot(yy[-which(is.na(avi.ES.s.full$Pol.u))],avi.ES.s$Pol.u)
plot(yy[-which(is.na(avi.ES.s.full$Pol.u))],avi.ES.s$Lin.u)
plot(yy2[-which(is.na(avi.ES.s.full$Pol.u))],avi.ES.s$Pol.u)
plot(yy2[-which(is.na(avi.ES.s.full$Pol.u))],avi.ES.s$Lin.u)

##################################################################
##################################################################
##################################################################


# function to understand the distribution of sample sizes for component studies

n.xy <- function (data) {
	
	xx <- split(data,list(data$Article_ID,data$Study_ID),drop=T)	
yy <- lapply(xx, function(y){
	dim(y)[1]})

yy2 <- do.call(rbind,yy)

return (yy2)
	
}

samp.sizes <- n.xy(data)

par(mfrow = c(1,2))
hist(samp.sizes,breaks = 100,xlim=c(0,500))
hist(rgamma(1115,shape = 0.3,scale = 100),xlim=c(0,500),breaks = 30)

median(samp.sizes)

write.table (samp.sizes, file = "AvI_Study_Sample_Sizes.txt", sep = "\t", quote = F, col.names = T, row.names = T, na = "")


##################################################################
##################################################################
##################################################################
## Plots of the raw data
#subset the data
xx <- data[which(data$Multi_spp_resp=="SINGLE"&data$Trophic_level=="Above"),]
#split the data by study
xx <- split(xx,xx$Art_stud, drop = T)
#plot the data
par(mfrow=c(4,5))
lapply(xx[61:76], function(y){
	plot(y$Abundance_Invader, y$Response)
})

#subset the data
xx <- data[which(data$Multi_spp_resp=="MULTIPLE"&data$Trophic_level=="Intra"),]
#split the data by study
xx <- split(xx,xx$Art_stud, drop = T)
# sub-sample to 100 randomly chosen studies
xx <- xx[sample(1:363, 100)]
#plot the data
par(mfrow=c(4,5))
lapply(xx[81:100], function(y){
	plot(y$Abundance_Invader, y$Response)
})

##################################################################
##################################################################
##################################################################