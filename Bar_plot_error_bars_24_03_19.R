# Code written by Raj Whitlock
# A series of functions to draw error whiskers/bars around points and bars on plots in R (both existing plots and plotting de novo), with control over formatting


# requires gtools

library(gtools)

##########################################################################
##########################################################################
##########################################################################

# plot a barplot and accompanying error bars
# this function is legacy code - I recommend users first draw their plots or barcharts, then add error bars to them afterwards using either e.bars.points2 or e.bars.points3

# xbar - bar height
# z - bar heights as offset from xbar
# barnames - bar names
# cap - relative width of cap of error bar
# ymin - logical, if true, selects a new x-axis intercept (minimum y value) based on the data
# ymax - logical, if true, selects a new y-axis maximum plotting extent based on the data
# plottitle - plot title
# xlab - same as par xlab
# ylab - same as par ylab
# col - vector of colours for bars
# txtpoints - controls text size via cex
# density, angle, width, space - as defined in the barplot function
# lty - same as par lty
# lower.bars - logical, if true then symmetrical lower bars are drawn

error.bars<-function(xbar,z,barnames,cap=0.3,ymin=F,ymax = "F", plottitle="",xlab="",ylab="",col="white",txtpoints=1,
density=NULL,angle=45,width=1,space=0.2,addline=F,lty=2,lower.bars=F,...){

if(ymin=="F"){
xdig<--(floor(log(min(xbar),10)))
ylower<-floor(10^xdig*min(xbar))/10^xdig}

if(ymin!="F"){
ylower<-ymin}

if(ymax == "F") {yupper <- max(xbar)+(1.5*max(z))} else {yupper <- ymax}
if((ylower<0)==T)ylower<-0
xv<-barplot(xbar-ylower,ylim=c(ylower,yupper),col=col,xpd=F,main=plottitle,xlab=xlab,ylab=ylab,
cex.axis=txtpoints,cex.lab=txtpoints,cex.names=txtpoints,tck=0.025,tcl=0.5,bty="l",offset=ylower,
names.arg=barnames,width=width,space=space,density=density,angle=angle, cex.main = 1.5, font.main = 3)
if(cap==0)g<-0
if(cap!=0)g<-width*cap
for(i in 1:length(xv)){
lines(c(xv[i],xv[i]),c(xbar[i]+z[i],xbar[i]))
lines(c(xv[i]-g,xv[i]+g),c(xbar[i]+z[i],xbar[i]+z[i]))
if(lower.bars==T){
lines(c(xv[i],xv[i]),c(xbar[i]-z[i],xbar[i]))
lines(c(xv[i]-g,xv[i]+g),c(xbar[i]-z[i],xbar[i]-z[i]))	
	}
if(is.numeric(addline)==T){
	lines(c(xv[i]-(width/1.7),xv[i]+(width/1.7)),c(addline[i],addline[i]),lty=lty,col="gray20")}

}}

##########################################################################
##########################################################################
##########################################################################

# plot a barplot and accompanying error bars, where upper and lower error bars can have asymmetrical length (simplified/ improved version of errorbars)
# this function is legacy code - I recommend users first draw their plots or barcharts, then add error bars to them afterwards using either e.bars.points2 or e.bars.points3

# yv - bar height
# z - upper error bar whisker lengths as offset from xbar
# xx - lower error bar whisker lengths as offset from xbar
# col - colours for error bars
# ylab - same as par ylab
# txtpoints - controls text size via cex 
# names.arg - bar names
# width - sets relative width of error bar cap, similar to cap, above
# ylim - same as par ylim
# offset - same as barplot offset
# ... other arguments that can be passed to barplot


error.bars2<-function(yv,z,xx = F,col=1,ylab="",txtpoints=1,names.arg="",width=10,ylim,offset=0,...){
xv<-barplot(yv,ylim=ylim,space=0.3,col=colors()[col],ylab=ylab,cex.axis=txtpoints,cex.lab=txtpoints,cex.names=txtpoints,tck=0.025,names.arg=names.arg,xpd=F,offset=offset)
g<-(max(xv)-min(xv))/width
if(xx[1]==F){
for(i in 1:length(xv)){
lines(c(xv[i],xv[i]),c(yv[i]+z[i],yv[i]))
lines(c(xv[i]-g,xv[i]+g),c(yv[i]+z[i],yv[i]+z[i]))
}}	else	{
	for(i in 1:length(xv)){
	lines(c(xv[i],xv[i]),c(yv[i]+z[i],yv[i]))
	lines(c(xv[i]-g,xv[i]+g),c(yv[i]+z[i],yv[i]+z[i]))
	lines(c(xv[i],xv[i]),c(yv[i]-xx[i],yv[i]))
	lines(c(xv[i]-g,xv[i]+g),c(yv[i]-xx[i],yv[i]-xx[i]))
	}}

}



##########################################################################
##########################################################################
##########################################################################

# plots simple HORIZONTAL ERROR BARS about existing points or bars on a plot

# bar - vector of bar heights
# x - vector of bar/ point positions on the y axis (centre of each bar/ point on the y axis)
# z1 - position of left error bar whisker, in absolute terms (not as an offset) on the x axis
# z2 - position of right error bar whisker, in absolute terms (not as an offset) on the x axis
# cap - relative width of cap/ whisker of error bar
# lty - same as par lty
# lwd - same as par lwd
# clr - vector of colours for error bars

e.bars.points2<-function(bar, x = NULL, z1,z2,cap=0.1,lwd=1,lty=1,clr){
bar <- as.numeric(bar)
x <- as.numeric(x)
z1 <- as.numeric(z1)
z2 <- as.numeric(z2)

		for(i in 1:length(bar)){

	lines(c(z1[i],bar[i]),c(x[i],x[i]),lwd=lwd,lty=lty,col=clr[i])
	lines(c(z1[i],z1[i]),c(x[i]-cap,x[i]+cap),lwd=lwd,lty=lty,col=clr[i])
	lines(c(z2[i],bar[i]),c(x[i],x[i]),lwd=lwd,lty=lty,col=clr[i])
	lines(c(z2[i],z2[i]),c(x[i]-cap,x[i]+cap),lwd=lwd,lty=lty,col=clr[i])
	}
	}

##########################################################################
##########################################################################
##########################################################################

# plots simple VERTICAL ERROR BARS about existing points or bars on a plot

# bar - vector of bar heights
# x - vector of bar/ point positions on the x axis (centre of each bar/ point on the x axis)
# z1 - position of upper error bar whisker, in absolute terms (not as an offset) on the y axis
# z2 - position of lower error bar whisker, in absolute terms (not as an offset) on the y axis
# cap - relative width of cap/ whisker of error bar
# lty - same as par lty
# lwd - same as par lwd
# clr - vector of colours for error bars
# lower - logical, if true, lower error bars are drawn

e.bars.points3<-function(bar, x = NULL, z1,z2=NA,cap=0.1,lwd=1,lty=1,clr,lower = F){

bar <- as.numeric(bar)
x <- as.numeric(x)
z1 <- as.numeric(z1)
z2 <- as.numeric(z2)
if(length(clr)==1)clr <- rep(clr,length(bar))

		for(i in 1:length(bar)){

	lines(c(x[i],x[i]),c(z1[i],bar[i]),lwd=lwd,lty=lty,col=clr[i])
	lines(c(x[i]-cap,x[i]+cap),c(z1[i],z1[i]),lwd=lwd,lty=lty,col=clr[i])
	if(lower==T) {lines(c(x[i],x[i]),c(z2[i],bar[i]),lwd=lwd,lty=lty,col=clr[i])
	lines(c(x[i]-cap,x[i]+cap),c(z2[i],z2[i]),lwd=lwd,lty=lty,col=clr[i])}
	}
	}


