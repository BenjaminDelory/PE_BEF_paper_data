############################################################################################################################
#R code used in Delory et al (2019) When history matters: the overlooked role of priority effects in grassland overyielding
############################################################################################################################

#########################################################################
#Set the path to the PE_BEF_paper_data directory
#Example: "C:/Users/Delory/Documents/GitHub/PE_BEF_paper_data"
#########################################################################

path<-"SET_WORKING_DIRECTORY_HERE"

##################
#Run the code :-)
##################

#Set working directory

setwd(path)

#Install bef R package from GitHub

install.packages("devtools")
devtools::install_github("BenjaminDelory/bef")

#Load R packages

library(bef)
library(readr)
library(beanplot)
library(emmeans)
library(smatr)
library(ggplot2)
library(ggpubr)
library(vegan)

#Create bootstrap function to calculate confidence intervals

CIbootstrap<-function(x, data=NULL, n=10000, CI=0.95){
  
  #x is a numeric vector containing sample values
  #x can also be a formula (y~factor1*factor2)
  #data is a data frame
  #n is the number of bootsrapping
  #CI is the confidence interval type
  
  if (class(x)=="formula"){
    
    response<-as.character(x[[2]])
    factor1<-as.character(x[[3]][[2]])
    factor2<-as.character(x[[3]][[3]])
    
    tmp<-expand.grid(factor1=levels(as.factor(data[,factor1])), factor2=levels(as.factor(data[,factor2])))
    colnames(tmp)<-c(factor1, factor2)
    tmp$mean<-NA
    tmp$CIlow<-NA
    tmp$CIhigh<-NA
    
    for (l in 1:nrow(tmp)){
      
      x<-data[which(data[,factor1]==tmp[l,factor1] & data[,factor2]==tmp[l,factor2]),response]
      if (TRUE %in% is.na(x)){x<-x[-which(is.na(x)==TRUE)]}
      bstrap<-c()
      
      for (i in 1:n){
        
        bsample<-sample(x, length(x), replace=TRUE)
        bestimate<-mean(bsample)
        bstrap<-c(bstrap, bestimate)}
      
      tmp[l,"mean"]<-mean(x)
      tmp[l,"CIlow"]<-quantile(bstrap, (1-CI)/2)
      tmp[l,"CIhigh"]<-quantile(bstrap, CI+(1-CI)/2)}
    
    return(tmp)}
  
  else {
  
    if (TRUE %in% is.na(x)){x<-x[-which(is.na(x)==TRUE)]}
    
    bstrap<-c()
    
    for (i in 1:n){
      
      bsample<-sample(x, length(x), replace=TRUE)
      bestimate<-mean(bsample)
      bstrap<-c(bstrap, bestimate)}
    
    return(c(quantile(bstrap, (1-CI)/2), quantile(bstrap, CI+(1-CI)/2)))}}

#Significance test based on 95CI

testH0<-function(x){
  
  x$test<-NA
  
  for (i in 1:nrow(x)){
    
    if (x$CIhigh[i]<x$H0[i]){x$test[i]<-"-"}
    if (x$CIlow[i]>x$H0[i]){x$test[i]<-"+"}
    if (x$CIlow[i]<x$H0[i] & x$CIhigh[i]>x$H0[i]){x$test[i]<-"0"}}
  
  return(x)}

#Calculate RYT (relative yield total)
RYT<-function(x){
  
  x$RY<-x$SDW/x$Mono
  results<-aggregate(x[,"RY"], by=list(Plot=x$Plot, Arrival=x$Arrival), sum)
  return(results)}

#Load data
data2013 <- read_delim("Data/Data_PE_Julich_June2013.txt",
                   "\t", escape_double = FALSE, trim_ws = TRUE)
View(data2013)
data2013<-as.data.frame(data2013)

data2013$SDW<-data2013$SDW*10

mono2013<-data2013[data2013$Plot=="Mono",]
mix2013<-data2013[data2013$Plot!="Mono",]

mix2013<-aggregate(mix2013$SDW, by=list(Art=mix2013$Art, Plot=mix2013$Plot), mean)

mix1<-matrix(ncol=9, nrow=32)
colnames(mix1)<-unique(mix2013$Art)
rownames(mix1)<-unique(mix2013$Plot)

for (i in 1:nrow(mix2013)){mix1[mix2013$Plot[i],mix2013$Art[i]]<-mix2013$x[i]}
mix1[is.na(mix1)==TRUE]<-0
mix1<-mix1[c(1:4,9:12,17:20,25:28),]
mix1.2013<-mix1

mono1<-matrix(mono2013$SDW, ncol=9, nrow=4, byrow=FALSE)
colnames(mono1)<-unique(mono2013$Art)
rownames(mono1)<-c("Mono1", "Mono2", "Mono3", "Mono4")
mono1.2013<-mono1

#Additive partitioning
res2013<-apm(mix=mix1, mono=mono1, method="fox")
res2013<-as.data.frame(res2013)
res2013$Arrival<-c(rep("F-first", 4), rep("G-first", 4), rep("aSynch", 4), rep("L-first", 4))
res2013$CE<-res2013$TDCE+res2013$TICE
res2013$Year<-rep(2013, nrow(res2013))
res2013$Plot<-rownames(res2013)
res2013$SE<-res2013$TDCE+res2013$DE

#Prepare data 2014
data2014 <- read_delim("Data/Data_PE_Julich_June2014.txt",
"\t", escape_double = FALSE, trim_ws = TRUE)
View(data2014)
data2014<-as.data.frame(data2014)

data2014<-data2014[-which(data2014$Art=="Arrela"|data2014$Art=="Poapra"|data2014$Art=="Trirep"),]

mono2014<-data2014[data2014$Plot=="Mono",]
mix2014<-data2014[data2014$Plot!="Mono",]

mix2014<-aggregate(mix2014$SDW, by=list(Art=mix2014$Art, Plot=mix2014$Plot), mean)

mix1<-matrix(ncol=9, nrow=32)
colnames(mix1)<-unique(mix2014$Art)
rownames(mix1)<-unique(mix2014$Plot)

for (i in 1:nrow(mix2014)){mix1[mix2014$Plot[i],mix2014$Art[i]]<-mix2014$x[i]}
mix1[is.na(mix1)==TRUE]<-0
mix1<-mix1[c(1:4,9:12,17:20,25:28),]
mix1.2014<-mix1

mono1<-matrix(mono2014$SDW, ncol=9, nrow=4, byrow=FALSE)
colnames(mono1)<-unique(mono2014$Art)
rownames(mono1)<-c("Mono1", "Mono2", "Mono3", "Mono4")
mono1.2014<-mono1

#Additive partitioning 2014
res2014<-apm(mix=mix1, mono=mono1, method="fox")
res2014<-as.data.frame(res2014)
res2014$Arrival<-c(rep("aSynch", 4), rep("F-first", 4), rep("G-first", 4), rep("L-first", 4))
res2014$CE<-res2014$TDCE+res2014$TICE
res2014$Year<-rep(2014, nrow(res2014))
res2014$Plot<-rownames(res2014)
res2014$Arrival<-as.factor(res2014$Arrival)
res2014$SE<-res2014$TDCE+res2014$DE

#Merge data 2014 and 2013
res<-rbind(res2013, res2014)
res$Year<-factor(res$Year)
res$Arrival<-factor(res$Arrival)
rownames(res)<-NULL

res$rTICE<-(100/res$NBE)*res$TICE
res$rTDCE<-(100/res$NBE)*res$TDCE
res$rDE<-(100/res$NBE)*res$DE

#NBE
model<-lm(NBE~Year*Arrival, data=res)
plot(fitted(model), resid(model), pch=16); abline(h=0, lty=2)
anova(model) #Only a year effect in net biodiversity effect
capture.output(anova(model), file="NBE.csv")

posthoc<-emmeans(model, "Year", contr="tukey", type="response") #Year effect
summary(posthoc, adjust="tukey") #Get groups
posthoc<-emmeans(model, "Arrival", by="Year", contr="tukey", type="response") #Analyse interaction
summary(posthoc, adjust="tukey") #Get groups

#DE
model<-lm(DE~Year*Arrival, data=res)
plot(fitted(model), resid(model), pch=16); abline(h=0, lty=2)
anova(model) #Significant Year*Arrival interaction
capture.output(anova(model), file="DE.csv")

posthoc<-emmeans(model, "Arrival", by="Year", contr="tukey", type="response") #Analyse interaction
summary(posthoc, adjust="tukey") #Get groups

#TICE
model<-lm(TICE~Year*Arrival, data=res)
plot(fitted(model), resid(model), pch=16); abline(h=0, lty=2)
anova(model) #Significant Year*Arrival interaction
capture.output(anova(model), file="TICE.csv")

posthoc<-emmeans(model, "Arrival", by="Year", contr="tukey", type="response") #Analyse interaction
summary(posthoc, adjust="tukey") #Get groups

#TDCE
model<-lm(TDCE~Year*Arrival, data=res)
plot(fitted(model), resid(model), pch=16); abline(h=0, lty=2)
anova(model) #Significant Year effect
capture.output(anova(model), file="TDCE.csv")

posthoc<-emmeans(model, "Year", contr="tukey", type="response") #Year effect
summary(posthoc, adjust="tukey") #Get groups
posthoc<-emmeans(model, "Arrival", by="Year", contr="tukey", type="response") #Analyse interaction
summary(posthoc, adjust="tukey") #Get groups

#CE
model<-lm(CE~Year*Arrival, data=res)
plot(fitted(model), resid(model), pch=16); abline(h=0, lty=2)
anova(model) #Significant interaction effect

posthoc<-emmeans(model, "Arrival", by="Year", contr="tukey", type="response") #Analyse interaction
summary(posthoc, adjust="tukey") #Get groups

#SE
model<-lm(SE~Year*Arrival, data=res)
plot(fitted(model), resid(model), pch=16); abline(h=0, lty=2)
anova(model)

posthoc<-emmeans(model, "Arrival", by="Year", contr="tukey", type="response") #Analyse interaction
summary(posthoc, adjust="tukey") #Get groups

###############
#Plot figure 2
###############

res<-as.data.frame(res)
res2013<-as.data.frame(res2013)
res2014<-as.data.frame(res2014)
colors<-c("white", "grey50", "black", "black")
x<-c(1,2.25,3.5,4.75,6.75,8,9.25,10.5)
res$x<-c(rep(x[2],4), rep(x[3],4), rep(x[1],4), rep(x[4],4),
         rep(x[5],4), rep(x[6],4), rep(x[7],4), rep(x[8],4))-0.15
x<-x+0.15

tiff(filename="Figure2.tif", res=600, height=15, width=16, units="cm", compression="lzw", pointsize=8)

#NBE
par(bty="l", mar=c(2,4.5,4,2)+0,1, fig=c(0,0.5,0.5,1))
gf<-beanplot(NBE~Arrival*Year, res, log="", what=c(0,0,0,0), border = NA, 
             ylab=expression(paste("Net biodiversity effect (g ", m^-2, ")", sep="")), xlab="", 
             las=1, ylim=c(-200,600), cex.lab=1, cex.axis=1, 
             names=rep(c("", "", "", ""), 2), at=x,
             xlim=c(0.5,11), main="", yaxt="n")
axis(2, at=seq(-200, 600, by=100), las=1)
points(res$x, res$NBE, pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(h=0, lty=2)
model<-lm(NBE~Year*Arrival, data=res)
posthoc<-emmeans(model, "Year", contr="tukey", type="response") #Year effect
summaryNBE<-summary(posthoc$emmeans, adjust="tukey")

segments(x0=0.5, x1=5.25, y0=summaryNBE$emmean[1], y1=summaryNBE$emmean[1], lty=3)
segments(x0=6.25, x1=11, y0=summaryNBE$emmean[2], y1=summaryNBE$emmean[2], lty=3)

posthoc<-emmeans(model, "Arrival", by="Year", contr="tukey", type="response") #Interaction effect
summaryNBE<-summary(posthoc$emmeans, adjust="tukey")
arrows(x0=x, y0=summaryNBE$lower.CL, x1=x, y1=summaryNBE$upper.CL, code=3, angle=90, length=0.02)
points(x, summaryNBE$emmean, pch=21, bg=c(rep("black",4), rep("white",3), "black"), col="black",
       cex=1.1)

text((2.25+3.5)/2, 600, "ns", cex=1)
text((8+9.25)/2, 350, "ns", cex=1)
text(0.6, 600, "a", cex=1.4, font=2)
text(11, 600, expression(paste(italic(P[Arrival]), " = 0.911")), pos=2, cex=0.9)
text(11, 550, expression(paste(italic(P[Year]), " < 0.001")), pos=2, cex=0.9)
text(11, 500, expression(paste(italic(P[Interaction]), " = 0.139")), pos=2, cex=0.9)

#DE
par(bty="l", mar=c(2,4.5,4,2)+0,1, fig=c(0.5,1,0.5,1), new=T)
gf<-beanplot(DE~Arrival*Year, res, log="", what=c(0,0,0,0), border = NA, 
             ylab=expression(paste("Dominance effect (g ", m^-2, ")", sep="")), xlab="", 
             las=1, ylim=c(-200,600), cex.lab=1, cex.axis=1, 
             names=rep(c("", "", "", ""), 2), at=x,
             xlim=c(0.5,11), main="", yaxt="n")
axis(2, at=seq(-200, 600, by=100), las=1)
points(res$x, res$DE, pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(h=0, lty=2)
model<-lm(DE~Year*Arrival, data=res)
posthoc<-emmeans(model, "Year", contr="tukey", type="response") #Year effect
summaryDE<-summary(posthoc$emmeans, adjust="tukey")

segments(x0=0.5, x1=5.25, y0=summaryDE$emmean[1], y1=summaryDE$emmean[1], lty=3)
segments(x0=6.25, x1=11, y0=summaryDE$emmean[2], y1=summaryDE$emmean[2], lty=3)

posthoc<-emmeans(model, "Arrival", by="Year", contr="tukey", type="response") #Interaction effect
summaryDE<-summary(posthoc$emmeans, adjust="tukey")
arrows(x0=x, y0=summaryDE$lower.CL, x1=x, y1=summaryDE$upper.CL, code=3, angle=90, length=0.02)
points(x, summaryDE$emmean, pch=21, bg=c(rep("black",4), rep("white",2), "black", "white"), 
       col="black", cex=1.1)

text(c(1,2.25,3.5,4.75)+0.15, summaryDE$upper.CL[1:4], c("a", "a", "a", "b"), cex=1, pos=3)
text((8+9.25)/2, 200, "ns", cex=1)
text(0.6, 600, "b", cex=1.4, font=2)
text(11, 600, expression(paste(italic(P[Arrival]), " < 0.001")), pos=2, cex=0.9)
text(11, 550, expression(paste(italic(P[Year]), " < 0.001")), pos=2, cex=0.9)
text(11, 500, expression(paste(italic(P[Interaction]), " = 0.008")), pos=2, cex=0.9)

#TICE
par(bty="l", mar=c(5,4.5,1,2)+0,1, fig=c(0,0.5,0,0.5), new=T)
gf<-beanplot(TICE~Arrival*Year, res, log="", what=c(0,0,0,0), border = NA, 
             ylab=expression(paste("Trait-independent complementarity effect (g ", m^-2, ")", sep="")), xlab="", 
             las=1, ylim=c(-200,600), cex.lab=1, cex.axis=1, 
             names=rep(c("Sync", expression(F[first]), expression(G[first]), expression(L[first])), 2), at=x,
             xlim=c(0.5,11), main="", yaxt="n")
axis(2, at=seq(-200, 600, by=100), las=1)
points(res$x, res$TICE, pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(h=0, lty=2)
model<-lm(TICE~Year*Arrival, data=res)
posthoc<-emmeans(model, "Year", contr="tukey", type="response") #Year effect
summaryTICE<-summary(posthoc$emmeans, adjust="tukey")

segments(x0=0.5, x1=5.25, y0=summaryTICE$emmean[1], y1=summaryTICE$emmean[1], lty=3)
segments(x0=6.25, x1=11, y0=summaryTICE$emmean[2], y1=summaryTICE$emmean[2], lty=3)

posthoc<-emmeans(model, "Arrival", by="Year", contr="tukey", type="response") #Interaction effect
summaryTICE<-summary(posthoc$emmeans, adjust="tukey")
arrows(x0=x, y0=summaryTICE$lower.CL, x1=x, y1=summaryTICE$upper.CL, code=3, angle=90, length=0.02)
points(x, summaryTICE$emmean, pch=21, bg=c("black", "white", "black", "white", rep("white",4)), 
       col="black", cex=1.1)

text(c(1,2.25,3.5,4.75)+0.15, summaryTICE$upper.CL[1:4], c("a", "ab", "a", "b"), cex=1, pos=3)
text((8+9.25)/2, 300, "ns", cex=1)
text(0.6, 600, "c", cex=1.4, font=2)

text((2.25+3.5)/2, par("usr")[[3]]-170, "2013", cex=1.1, font=2, xpd=TRUE)
text((8+9.25)/2, par("usr")[[3]]-170, "2014", cex=1.1, font=2, xpd=TRUE)
text(11, 600, expression(paste(italic(P[Arrival]), " = 0.411")), pos=2, cex=0.9)
text(11, 550, expression(paste(italic(P[Year]), " = 0.021")), pos=2, cex=0.9)
text(11, 500, expression(paste(italic(P[Interaction]), " = 0.011")), pos=2, cex=0.9)

#TDCE
par(bty="l", mar=c(5,4.5,1,2)+0,1, fig=c(0.5,1,0,0.5), new=T)
gf<-beanplot(TDCE~Arrival*Year, res, log="", what=c(0,0,0,0), border = NA, 
             ylab=expression(paste("Trait-dependent complementarity effect (g ", m^-2, ")", sep="")), xlab="", 
             las=1, ylim=c(-200,600), cex.lab=1, cex.axis=1, 
             names=rep(c("Sync", expression(F[first]), expression(G[first]), expression(L[first])), 2), at=x,
             xlim=c(0.5,11), main="", yaxt="n")
axis(2, at=seq(-200, 600, by=100), las=1)
points(res$x, res$TDCE, pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(h=0, lty=2)
model<-lm(TDCE~Year*Arrival, data=res)
posthoc<-emmeans(model, "Year", contr="tukey", type="response") #Year effect
summaryTDCE<-summary(posthoc$emmeans, adjust="tukey")

segments(x0=0.5, x1=5.25, y0=summaryTDCE$emmean[1], y1=summaryTDCE$emmean[1], lty=3)
segments(x0=6.25, x1=11, y0=summaryTDCE$emmean[2], y1=summaryTDCE$emmean[2], lty=3)

posthoc<-emmeans(model, "Arrival", by="Year", contr="tukey", type="response") #Interaction effect
summaryTDCE<-summary(posthoc$emmeans, adjust="tukey")
arrows(x0=x, y0=summaryTDCE$lower.CL, x1=x, y1=summaryTDCE$upper.CL, code=3, angle=90, length=0.02)
points(x, summaryTDCE$emmean, pch=21, bg=c("white", "white", "black", "white", rep("white",4)), 
       col="black", cex=1.1)

text((2.25+3.5)/2, 150, "ns", cex=1)
text((8+9.25)/2, 150, "ns", cex=1)
text(0.6, 600, "d", cex=1.4, font=2)

text((2.25+3.5)/2, par("usr")[[3]]-170, "2013", cex=1.1, font=2, xpd=TRUE)
text((8+9.25)/2, par("usr")[[3]]-170, "2014", cex=1.1, font=2, xpd=TRUE)
text(11, 600, expression(paste(italic(P[Arrival]), " = 0.698")), pos=2, cex=0.9)
text(11, 550, expression(paste(italic(P[Year]), " = 0.008")), pos=2, cex=0.9)
text(11, 500, expression(paste(italic(P[Interaction]), " = 0.070")), pos=2, cex=0.9)

dev.off()

#############################
#Plot supplementary figure 2
#############################

mono2013<-data2013[data2013$Plot=="Mono",]
mix2013<-data2013[data2013$Plot!="Mono",]

mix2013<-aggregate(mix2013$SDW, by=list(Art=mix2013$Art, Plot=mix2013$Plot), mean)

mix1<-matrix(ncol=9, nrow=32)
colnames(mix1)<-unique(mix2013$Art)
rownames(mix1)<-unique(mix2013$Plot)

for (i in 1:nrow(mix2013)){mix1[mix2013$Plot[i],mix2013$Art[i]]<-mix2013$x[i]}
mix1[is.na(mix1)==TRUE]<-0
mix1<-mix1[c(1:4,9:12,17:20,25:28),]

sp2013<-mix1
sp2013<-sp2013[order(rownames(sp2013)),]
sp2013<-sp2013[,order(colnames(sp2013))]
mix2013<-data.frame(SDW=as.vector(sp2013),
                    Art=c(rep("Achmil",16),rep("Dacglo",16),rep("Fespra",16),rep("Hollan",16),rep("Leuvul",16),
                          rep("Lotcor",16),rep("Medsat",16),rep("Plalan",16),rep("Tripra",16)),
                    Arrival=rep(c(rep("F-first",4),rep("G-first",4),rep("Sync",4),rep("L-first",4)), 9),
                    Plot=rep(rownames(sp2013), 9))

mono2013<-data2013[data2013$Plot=="Mono",]
mono2013<-aggregate(mono2013$SDW, by=list(Art=mono2013$Art), mean)
mix2013$Mono<-mono2013$x[match(mix2013$Art, mono2013$Art)]
mono2013$pch<-c(21,22,22,22,21,24,24,21,24)
pal<-c("#009E73", "#F0E442", "#0072B2")
mono2013$bg<-c(pal[3],pal[2],pal[2],pal[2],pal[3],pal[1],pal[1],pal[3],pal[1])

RYT2013<-RYT(mix2013)
beanplot(x~Arrival, RYT2013)

CI<-CIbootstrap(SDW~Art*Arrival, mix2013)
CI$Mono<-mono2013$x[match(CI$Art, mono2013$Art)]
CI$H0<-CI$Mono/9
testH0(CI)

tiff(filename="FigureS2.tif", res=600, compression="lzw", width=15, height=13, units="cm", 
     pointsize=7)

par(mar=c(2.5,4.5,4,1)+0.1, fig=c(0,0.5,0.5,1))
plot(mix2013$Mono[mix2013$Arrival=="Sync"], mix2013$SDW[mix2013$Arrival=="Sync"], 
     xlab="", 
     ylab=expression(paste("Yield in mixture (g ", m^-2, ")", sep="")), las=1, type="n",
     xaxt="n", yaxt="n", xlim=c(150,750), ylim=c(0,800))
axis(1, at=seq(from=150, to=750, by=100), las=1)
axis(2, at=seq(0,800, by=100), las=1)
axis(3, at=mono2013$x[-c(3,6,9)], labels=mono2013$Art[-c(3,6,9)], las=2, cex.axis=0.8, padj=c(0,1,rep(0.5,4)))
axis(3, at=mono2013$x[c(3,6,9)], labels=mono2013$Art[c(3,6,9)], las=2, cex.axis=0.8, padj=rep(0.5,3), font=2)
points(mix2013$Mono[mix2013$Arrival=="Sync"], mix2013$SDW[mix2013$Arrival=="Sync"], pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(a=0, b=1, lty=1)
lines(c(1:1000), c(1:1000)/9, lty=2)
lines(c(mean(mono2013$x), mean(mono2013$x)), c(-100,mean(mono2013$x)), lty=3)
arrows(x0=mono2013$x, y0=CI$CIlow[CI$Arrival=="Sync"], 
       x1=mono2013$x, y1=CI$CIhigh[CI$Arrival=="Sync"],
       code=3, length=0.02, angle=90)
points(mono2013$x, CI$mean[CI$Arrival=="Sync"], pch=mono2013$pch, col="black", 
       bg=mono2013$bg, cex=1.2)
text(160,770, "Synchronous", cex=1, pos=4)
text(120,770, "a", cex=1.4, font=2, pos=4)
legend(130,720, legend=c("Mixture=Monoculture", "Null model", "Average monoculture yield"), lty=c(1,2,3), bty="n", cex=0.8)

par(mar=c(2.5,3.5,4,2)+0.1, fig=c(0.5,1,0.5,1), new=TRUE)
plot(mix2013$Mono[mix2013$Arrival=="F-first"], mix2013$SDW[mix2013$Arrival=="F-first"], 
     xlab="", 
     ylab="", las=1, type="n",
     xaxt="n", yaxt="n", xlim=c(150,750), ylim=c(0,800))
axis(1, at=seq(from=150, to=750, by=100), las=1)
axis(2, at=seq(0,800, by=100), las=1)
axis(3, at=mono2013$x[-c(3,6,9)], labels=mono2013$Art[-c(3,6,9)], las=2, cex.axis=0.8, padj=c(0,1,rep(0.5,4)))
axis(3, at=mono2013$x[c(3,6,9)], labels=mono2013$Art[c(3,6,9)], las=2, cex.axis=0.8, padj=rep(0.5,3), font=2)
points(mix2013$Mono[mix2013$Arrival=="F-first"], mix2013$SDW[mix2013$Arrival=="F-first"], pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(a=0, b=1, lty=1)
lines(c(1:1000), c(1:1000)/9, lty=2)
lines(c(mean(mono2013$x), mean(mono2013$x)), c(-100,mean(mono2013$x)), lty=3)
arrows(x0=mono2013$x, y0=CI$CIlow[CI$Arrival=="F-first"], 
       x1=mono2013$x, y1=CI$CIhigh[CI$Arrival=="F-first"],
       code=3, length=0.02, angle=90)
points(mono2013$x, CI$mean[CI$Arrival=="F-first"], pch=mono2013$pch, col="black", 
       bg=mono2013$bg, cex=1.2)
text(160,770, "F-first", cex=1, pos=4)
text(120,770, "b", cex=1.4, font=2, pos=4)
legend(130,720, legend=c("Forbs", "Grasses", "Legumes"), pch=c(21,22,24), col="black", pt.bg=c(pal[3], pal[2], pal[1]), cex=0.9,
       bty="n", pt.cex=1.2)

par(mar=c(5,4.5,1.5,1)+0.1, fig=c(0,0.5,0,0.5), new=TRUE)
plot(mix2013$Mono[mix2013$Arrival=="G-first"], mix2013$SDW[mix2013$Arrival=="G-first"], 
     xlab=expression(paste("Yield in monoculture (g ", m^-2, ")", sep="")), 
     ylab=expression(paste("Yield in mixture (g ", m^-2, ")", sep="")), las=1, type="n",
     xaxt="n", yaxt="n", xlim=c(150,750), ylim=c(0,800))
axis(1, at=seq(from=150, to=750, by=100), las=1)
axis(2, at=seq(0,800, by=100), las=1)
axis(3, at=mono2013$x, labels=FALSE, las=2, cex.axis=0.8, padj=c(0,1,rep(0.5,7)))
points(mix2013$Mono[mix2013$Arrival=="G-first"], mix2013$SDW[mix2013$Arrival=="G-first"], pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(a=0, b=1, lty=1)
lines(c(mean(mono2013$x), mean(mono2013$x)), c(-100,mean(mono2013$x)), lty=3)
lines(c(1:1000), c(1:1000)/9, lty=2)
arrows(x0=mono2013$x, y0=CI$CIlow[CI$Arrival=="G-first"], 
       x1=mono2013$x, y1=CI$CIhigh[CI$Arrival=="G-first"],
       code=3, length=0.02, angle=90)
points(mono2013$x, CI$mean[CI$Arrival=="G-first"], pch=mono2013$pch, col="black", 
       bg=mono2013$bg, cex=1.2)
text(160,770, "G-first", cex=1, pos=4)
text(120,770, "c", cex=1.4, font=2, pos=4)

par(mar=c(5,3.5,1.5,2)+0.1, fig=c(0.5,1,0,0.5), new=TRUE)
plot(mix2013$Mono[mix2013$Arrival=="L-first"], mix2013$SDW[mix2013$Arrival=="L-first"], 
     xlab=expression(paste("Yield in monoculture (g ", m^-2, ")", sep="")), 
     ylab="", las=1, type="n",
     xaxt="n", yaxt="n", xlim=c(150,750), ylim=c(0,800))
axis(1, at=seq(from=150, to=750, by=100), las=1)
axis(2, at=seq(0,800, by=100), las=1)
axis(3, at=mono2013$x, labels=FALSE, las=2, cex.axis=0.8, padj=c(0,1,rep(0.5,7)))
points(mix2013$Mono[mix2013$Arrival=="L-first"], mix2013$SDW[mix2013$Arrival=="L-first"], pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(a=0, b=1, lty=1)
lines(c(1:1000), c(1:1000)/9, lty=2)
lines(c(mean(mono2013$x), mean(mono2013$x)), c(-100,mean(mono2013$x)), lty=3)
arrows(x0=mono2013$x, y0=CI$CIlow[CI$Arrival=="L-first"], 
       x1=mono2013$x, y1=CI$CIhigh[CI$Arrival=="L-first"],
       code=3, length=0.02, angle=90)
points(mono2013$x, CI$mean[CI$Arrival=="L-first"], pch=mono2013$pch, col="black", 
       bg=mono2013$bg, cex=1.2)
text(160,770, "L-first", cex=1, pos=4)
text(120,770, "d", cex=1.4, font=2, pos=4)

dev.off()

#############################
#Plot supplementary figure 3
#############################

mono2014<-data2014[data2014$Plot=="Mono",]
mix2014<-data2014[data2014$Plot!="Mono",]

mix2014<-aggregate(mix2014$SDW, by=list(Art=mix2014$Art, Plot=mix2014$Plot), mean)

mix1<-matrix(ncol=9, nrow=32)
colnames(mix1)<-unique(mix2014$Art)
rownames(mix1)<-unique(mix2014$Plot)

for (i in 1:nrow(mix2014)){mix1[mix2014$Plot[i],mix2014$Art[i]]<-mix2014$x[i]}
mix1[is.na(mix1)==TRUE]<-0
mix1<-mix1[c(1:4,9:12,17:20,25:28),]

sp2014<-mix1
sp2014<-sp2014[order(rownames(sp2014)),]
sp2014<-sp2014[,order(colnames(sp2014))]
mix2014<-data.frame(SDW=as.vector(sp2014),
                    Art=c(rep("Achmil",16),rep("Dacglo",16),rep("Fespra",16),rep("Hollan",16),rep("Leuvul",16),
                          rep("Lotcor",16),rep("Medsat",16),rep("Plalan",16),rep("Tripra",16)),
                    Arrival=rep(c(rep("Sync",4),rep("F-first",4),rep("G-first",4),rep("L-first",4)), 9),
                    Plot=rep(rownames(sp2014), 9))

mono2014<-data2014[data2014$Plot=="Mono",]
mono2014<-aggregate(mono2014$SDW, by=list(Art=mono2014$Art), mean)
mix2014$Mono<-mono2014$x[match(mix2014$Art, mono2014$Art)]
mono2014$pch<-c(21,22,22,22,21,24,24,21,24)
pal<-c("#009E73", "#F0E442", "#0072B2")
mono2014$bg<-c(pal[3],pal[2],pal[2],pal[2],pal[3],pal[1],pal[1],pal[3],pal[1])

RYT2014<-RYT(mix2014)
beanplot(x~Arrival, RYT2014)

CI<-CIbootstrap(SDW~Art*Arrival, mix2014)
CI$Mono<-mono2014$x[match(CI$Art, mono2014$Art)]
CI$H0<-CI$Mono/9
testH0(CI)

tiff(filename="FigureS3.tif", res=600, compression="lzw", width=15, height=13, units="cm", 
     pointsize=7)

par(mar=c(2.5,4.5,4,1)+0.1, fig=c(0,0.5,0.5,1))
plot(mix2014$Mono[mix2014$Arrival=="Sync"], mix2014$SDW[mix2014$Arrival=="Sync"], 
     xlab="", 
     ylab=expression(paste("Yield in mixture (g ", m^-2, ")", sep="")), las=1, type="n",
     xaxt="n", yaxt="n", xlim=c(50,500), ylim=c(0,500))
axis(1, at=seq(from=50, to=500, by=100), las=1)
axis(2, at=seq(0,500, by=100), las=1)
axis(3, at=mono2014$x[-c(3:5)], labels=mono2014$Art[-c(3:5)], las=2, cex.axis=0.8, padj=c(0,1,0.5,0.5,0,1))
axis(3, at=mono2014$x[c(3:5)], labels=mono2014$Art[c(3:5)], las=2, cex.axis=0.8, padj=c(0,1,0.5), font=2)
points(mix2014$Mono[mix2014$Arrival=="Sync"], mix2014$SDW[mix2014$Arrival=="Sync"], pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(a=0, b=1, lty=1)
lines(c(1:1000), c(1:1000)/9, lty=2)
lines(c(mean(mono2014$x), mean(mono2014$x)), c(-100,mean(mono2014$x)), lty=3)
arrows(x0=mono2014$x, y0=CI$CIlow[CI$Arrival=="Sync"], 
       x1=mono2014$x, y1=CI$CIhigh[CI$Arrival=="Sync"],
       code=3, length=0.02, angle=90)
points(mono2014$x, CI$mean[CI$Arrival=="Sync"], pch=mono2014$pch, col="black", 
       bg=mono2014$bg, cex=1.2)
text(60,490, "Synchronous", cex=1, pos=4)
text(30,490, "a", cex=1.4, font=2, pos=4)
legend(35,465, legend=c("Mixture=Monoculture", "Null model", "Average monoculture yield"), lty=c(1,2,3), bty="n", cex=0.8)

par(mar=c(2.5,3.5,4,2)+0.1, fig=c(0.5,1,0.5,1), new=TRUE)
plot(mix2014$Mono[mix2014$Arrival=="F-first"], mix2014$SDW[mix2014$Arrival=="F-first"], 
     xlab="", 
     ylab="", las=1, type="n",
     xaxt="n", yaxt="n", xlim=c(50,500), ylim=c(0,500))
axis(1, at=seq(from=50, to=500, by=100), las=1)
axis(2, at=seq(0,500, by=100), las=1)
axis(3, at=mono2014$x[-c(3:5)], labels=mono2014$Art[-c(3:5)], las=2, cex.axis=0.8, padj=c(0,1,0.5,0.5,0,1))
axis(3, at=mono2014$x[c(3:5)], labels=mono2014$Art[c(3:5)], las=2, cex.axis=0.8, padj=c(0,1,0.5), font=2)
points(mix2014$Mono[mix2014$Arrival=="F-first"], mix2014$SDW[mix2014$Arrival=="F-first"], pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(a=0, b=1, lty=1)
lines(c(1:1000), c(1:1000)/9, lty=2)
lines(c(mean(mono2014$x), mean(mono2014$x)), c(-100,mean(mono2014$x)), lty=3)
arrows(x0=mono2014$x, y0=CI$CIlow[CI$Arrival=="F-first"], 
       x1=mono2014$x, y1=CI$CIhigh[CI$Arrival=="F-first"],
       code=3, length=0.02, angle=90)
points(mono2014$x, CI$mean[CI$Arrival=="F-first"], pch=mono2014$pch, col="black", 
       bg=mono2014$bg, cex=1.2)
text(60,490, "F-first", cex=1, pos=4)
text(30,490, "b", cex=1.4, font=2, pos=4)
legend(35,465, legend=c("Forbs", "Grasses", "Legumes"), pch=c(21,22,24), col="black", pt.bg=c(pal[3], pal[2], pal[1]), cex=0.9,
       bty="n", pt.cex=1.2)

par(mar=c(5,4.5,1.5,1)+0.1, fig=c(0,0.5,0,0.5), new=TRUE)
plot(mix2014$Mono[mix2014$Arrival=="G-first"], mix2014$SDW[mix2014$Arrival=="G-first"], 
     xlab=expression(paste("Yield in monoculture (g ", m^-2, ")", sep="")), 
     ylab=expression(paste("Yield in mixture (g ", m^-2, ")", sep="")), las=1, type="n",
     xaxt="n", yaxt="n", xlim=c(50,500), ylim=c(0,500))
axis(1, at=seq(from=50, to=500, by=100), las=1)
axis(2, at=seq(0,500, by=100), las=1)
axis(3, at=mono2014$x, labels=FALSE, las=2, cex.axis=0.8, padj=c(0,1,rep(0.5,7)))
points(mix2014$Mono[mix2014$Arrival=="G-first"], mix2014$SDW[mix2014$Arrival=="G-first"], pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(a=0, b=1, lty=1)
lines(c(mean(mono2014$x), mean(mono2014$x)), c(-100,mean(mono2014$x)), lty=3)
lines(c(1:1000), c(1:1000)/9, lty=2)
arrows(x0=mono2014$x, y0=CI$CIlow[CI$Arrival=="G-first"], 
       x1=mono2014$x, y1=CI$CIhigh[CI$Arrival=="G-first"],
       code=3, length=0.02, angle=90)
points(mono2014$x, CI$mean[CI$Arrival=="G-first"], pch=mono2014$pch, col="black", 
       bg=mono2014$bg, cex=1.2)
text(60,490, "G-first", cex=1, pos=4)
text(30,490, "c", cex=1.4, font=2, pos=4)

par(mar=c(5,3.5,1.5,2)+0.1, fig=c(0.5,1,0,0.5), new=TRUE)
plot(mix2014$Mono[mix2014$Arrival=="L-first"], mix2014$SDW[mix2014$Arrival=="L-first"], 
     xlab=expression(paste("Yield in monoculture (g ", m^-2, ")", sep="")), 
     ylab="", las=1, type="n",
     xaxt="n", yaxt="n", xlim=c(50,500), ylim=c(0,500))
axis(1, at=seq(from=50, to=500, by=100), las=1)
axis(2, at=seq(0,500, by=100), las=1)
axis(3, at=mono2014$x, labels=FALSE, las=2, cex.axis=0.8, padj=c(0,1,rep(0.5,7)))
points(mix2014$Mono[mix2014$Arrival=="L-first"], mix2014$SDW[mix2014$Arrival=="L-first"], pch=16, col=adjustcolor("black", alpha.f=0.3), cex=0.8)
abline(a=0, b=1, lty=1)
lines(c(1:1000), c(1:1000)/9, lty=2)
lines(c(mean(mono2014$x), mean(mono2014$x)), c(-100,mean(mono2014$x)), lty=3)
arrows(x0=mono2014$x, y0=CI$CIlow[CI$Arrival=="L-first"], 
       x1=mono2014$x, y1=CI$CIhigh[CI$Arrival=="L-first"],
       code=3, length=0.02, angle=90)
points(mono2014$x, CI$mean[CI$Arrival=="L-first"], pch=mono2014$pch, col="black", 
       bg=mono2014$bg, cex=1.2)
text(60,490, "L-first", cex=1, pos=4)
text(30,490, "d", cex=1.4, font=2, pos=4)

dev.off()

###############
#Plot figure 3
###############

#Calculate P for 2013

mix2013<-data2013[data2013$Plot!="Mono",]
mix2013<-aggregate(mix2013$SDW, by=list(Art=mix2013$Art, Plot=mix2013$Plot), mean)
mix1<-matrix(ncol=9, nrow=32)
colnames(mix1)<-unique(mix2013$Art)
rownames(mix1)<-unique(mix2013$Plot)

for (i in 1:nrow(mix2013)){mix1[mix2013$Plot[i],mix2013$Art[i]]<-mix2013$x[i]}
mix1[is.na(mix1)==TRUE]<-0
mix1<-mix1[c(1:4,9:12,17:20,25:28),]

simpson<-diversity(x=mix1, index="simpson")

mono2013<-data2013[data2013$Plot=="Mono",]
mono1<-matrix(mono2013$SDW, ncol=9, nrow=4, byrow=FALSE)
colnames(mono1)<-unique(mono2013$Art)
rownames(mono1)<-c("Mono1", "Mono2", "Mono3", "Mono4")

res2013<-apm(mix=mix1, mono=mono1, method="fox")
res2013<-as.data.frame(res2013)
res2013$Arrival<-c(rep("F-first", 4), rep("G-first", 4), rep("aSynch", 4), rep("L-first", 4))
res2013$CE<-res2013$TDCE+res2013$TICE
res2013$Year<-rep(2013, nrow(res2013))
res2013$Plot<-rownames(res2013)

mix1<-as.data.frame(mix1)
mix1$Forbs<-mix1$Achmil+mix1$Leuvul+mix1$Plalan
mix1$Grasses<-mix1$Dacglo+mix1$Fespra+mix1$Hollan
mix1$Legumes<-mix1$Medsat+mix1$Lotcor+mix1$Tripra
mix1$Arrival<-c(rep("F-first",4), rep("G-first",4), rep("Sync",4), rep("L-first",4))
mix1$Replicate<-rep(1:4, 4)

pal<-c("white", "white", "white")
L<-res2013
L$Replicate<-rep(1:4, 4)
L$Simpson<-simpson
L$P<-NA
L$CIlow<-NA
L$CIhigh<-NA

for (i in 1:nrow(L)){
  
  plot<-rownames(L)[i]
  
  if (L$Arrival[i]=="F-first") {
    Ylast<-mix1$Legumes[which(rownames(mix1)==plot)]+mix1$Grasses[which(rownames(mix1)==plot)]
    Ysync1<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==1]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==1]
    Ysync2<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==2]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==2]
    Ysync3<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==3]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==3]
    Ysync4<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==4]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==4]}
  
  if (L$Arrival[i]=="G-first") {
    Ylast<-mix1$Legumes[which(rownames(mix1)==plot)]+mix1$Forbs[which(rownames(mix1)==plot)]
    Ysync1<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==1]+mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==1]
    Ysync2<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==2]+mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==2]
    Ysync3<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==3]+mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==3]
    Ysync4<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==4]+mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==4]}
  
  if (L$Arrival[i]=="L-first") {
    Ylast<-mix1$Forbs[which(rownames(mix1)==plot)]+mix1$Grasses[which(rownames(mix1)==plot)]
    Ysync1<-mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==1]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==1]
    Ysync2<-mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==2]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==2]
    Ysync3<-mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==3]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==3]
    Ysync4<-mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==4]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==4]}
  
  if (L$Arrival[i]!="aSynch"){
    
    P1<-2*(Ylast-Ysync1)/(Ysync1+abs(Ylast-Ysync1))
    P2<-2*(Ylast-Ysync2)/(Ysync2+abs(Ylast-Ysync2))
    P3<-2*(Ylast-Ysync3)/(Ysync3+abs(Ylast-Ysync3))
    P4<-2*(Ylast-Ysync4)/(Ysync4+abs(Ylast-Ysync4))
    
    L$P[i]<-(P1+P2+P3+P4)/4
    
    CI<-CIbootstrap(x=c(P1,P2,P3,P4), n=10000, CI=0.95)
    L$CIlow[i]<-CI[1]
    L$CIhigh[i]<-CI[2]}}

L<-L[-which(is.na(L$P)==TRUE),]
L2013<-L
P2013<-L$P

#Calculate P for 2014

mix2014<-data2014[data2014$Plot!="Mono",]
mix2014<-aggregate(mix2014$SDW, by=list(Art=mix2014$Art, Plot=mix2014$Plot), mean)
mix1<-matrix(ncol=9, nrow=32)
colnames(mix1)<-unique(mix2014$Art)
rownames(mix1)<-unique(mix2014$Plot)

for (i in 1:nrow(mix2014)){mix1[mix2014$Plot[i],mix2014$Art[i]]<-mix2014$x[i]}
mix1[is.na(mix1)==TRUE]<-0
mix1<-mix1[c(1:4,9:12,17:20,25:28),]

simpson<-diversity(x=mix1, index="simpson")

mono2014<-data2014[data2014$Plot=="Mono",]
mono1<-matrix(mono2014$SDW, ncol=9, nrow=4, byrow=FALSE)
colnames(mono1)<-unique(mono2014$Art)
rownames(mono1)<-c("Mono1", "Mono2", "Mono3", "Mono4")

res2014<-apm(mix=mix1, mono=mono1, method="fox")
res2014<-as.data.frame(res2014)
res2014$Arrival<-c(rep("aSynch", 4), rep("F-first", 4), rep("G-first", 4), rep("L-first", 4))
res2014$CE<-res2014$TDCE+res2014$TICE
res2014$Year<-rep(2014, nrow(res2014))
res2014$Plot<-rownames(res2014)

mix1<-as.data.frame(mix1)
mix1$Forbs<-mix1$Achmil+mix1$Leuvul+mix1$Plalan
mix1$Grasses<-mix1$Dacglo+mix1$Fespra+mix1$Hollan
mix1$Legumes<-mix1$Medsat+mix1$Lotcor+mix1$Tripra
mix1$Arrival<-c(rep("Sync", 4), rep("F-first", 4), rep("G-first", 4), rep("L-first", 4))
mix1$Replicate<-rep(1:4, 4)

pal<-c("white", "white", "white")
L<-res2014
L$Replicate<-rep(1:4, 4)
L$Simpson<-simpson
L$P<-NA
L$CIlow<-NA
L$CIhigh<-NA

for (i in 1:nrow(L)){
  
  plot<-rownames(L)[i]
  
  if (L$Arrival[i]=="F-first") {
    Ylast<-mix1$Legumes[which(rownames(mix1)==plot)]+mix1$Grasses[which(rownames(mix1)==plot)]
    Ysync1<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==1]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==1]
    Ysync2<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==2]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==2]
    Ysync3<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==3]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==3]
    Ysync4<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==4]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==4]}
  
  if (L$Arrival[i]=="G-first") {
    Ylast<-mix1$Legumes[which(rownames(mix1)==plot)]+mix1$Forbs[which(rownames(mix1)==plot)]
    Ysync1<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==1]+mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==1]
    Ysync2<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==2]+mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==2]
    Ysync3<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==3]+mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==3]
    Ysync4<-mix1$Legumes[mix1$Arrival=="Sync" & mix1$Replicate==4]+mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==4]}
  
  if (L$Arrival[i]=="L-first") {
    Ylast<-mix1$Forbs[which(rownames(mix1)==plot)]+mix1$Grasses[which(rownames(mix1)==plot)]
    Ysync1<-mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==1]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==1]
    Ysync2<-mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==2]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==2]
    Ysync3<-mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==3]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==3]
    Ysync4<-mix1$Forbs[mix1$Arrival=="Sync" & mix1$Replicate==4]+mix1$Grasses[mix1$Arrival=="Sync" & mix1$Replicate==4]}
  
  if (L$Arrival[i]!="aSynch"){
    
    P1<-2*(Ylast-Ysync1)/(Ysync1+abs(Ylast-Ysync1))
    P2<-2*(Ylast-Ysync2)/(Ysync2+abs(Ylast-Ysync2))
    P3<-2*(Ylast-Ysync3)/(Ysync3+abs(Ylast-Ysync3))
    P4<-2*(Ylast-Ysync4)/(Ysync4+abs(Ylast-Ysync4))
    
    L$P[i]<-(P1+P2+P3+P4)/4
    
    CI<-CIbootstrap(x=c(P1,P2,P3,P4), n=10000, CI=0.95)
    L$CIlow[i]<-CI[1]
    L$CIhigh[i]<-CI[2]}}

L<-L[-which(is.na(L$P)==TRUE),]
L2014<-L
P2014<-L$P

table<-rbind(L2013, L2014)
table$Year<-factor(table$Year)

#Relationship between P and Simpson diversity index
ggplot(table, aes(x=Simpson, y=P, color=Arrival))+
  geom_point(size=3)+
  facet_wrap(~Year)+
  theme_bw()+
  ylab(expression(paste("Priority effect index (", italic(P), ")", sep="")))+
  xlab("Simpson diversity index")

#Relationship between NBE and Simpson diversity index
ggplot(table, aes(x=Simpson, y=NBE, color=Arrival))+
  geom_point(size=3)+
  facet_wrap(~Year)+
  theme_bw()+
  ylab(expression(paste("Net biodiversity effect (g ", m^-2, ")")))+
  xlab("Simpson diversity index")

#Fit SMA models

modelNBE<-sma(NBE~P*Year, data=table, log="")
summary(modelNBE)
modelCE<-sma(CE~P*Year, data=table, log="")
summary(modelCE)
modelTICE<-sma(TICE~P*Year, data=table, log="")
summary(modelTICE)
modelTDCE<-sma(TDCE~P*Year, data=table, log="")
summary(modelTDCE)
modelDE<-sma(DE~P*Year, data=table, log="")
summary(modelDE)

#Plot figure 3

cbPalette <- c("#0072B2", "#F0E442", "#009E73")

cex.axis=10
cex.title=10

model<-data.frame(
  newx=seq(min(L2013$P),max(L2013$P),by=0.01),
  prd=modelNBE$coef$`2013`$`coef(SMA)`[2]*seq(min(L2013$P),max(L2013$P),by=0.01)+modelNBE$coef$`2013`$`coef(SMA)`[1])

p1<-ggplot(L2013, aes(x=P, y=NBE, fill=Arrival, shape=Arrival))+
  geom_segment(aes(x=CIlow, y=NBE, xend=CIhigh, yend=NBE), color="gray50")+
  geom_point(color="black", size=2)+
  geom_line(aes(x=newx, y=prd), data=model, color="black", size=1, inherit.aes = FALSE)+
  theme_bw()+
  facet_wrap(~Year)+
  scale_y_continuous(breaks=seq(0,600, by=100), limits=c(0,600))+
  geom_vline(xintercept = 0, linetype=3, color="black")+
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.8, units="lines"),
        legend.position=c(0.85,0.2),
        legend.background = element_rect(fill=NA),
        panel.grid=element_blank(),
        axis.text.x=element_text(size=cex.axis,color="black"),
        axis.text.y=element_text(size=cex.axis,color="black"),
        axis.title.y=element_text(size=cex.title, margin=margin(t=0,r=7,b=0,l=0)),
        axis.title.x=element_text(size=cex.title, margin=margin(t=7,r=0,b=0,l=0)),
        strip.text=element_text(size=12),
        panel.border = element_rect(color="black"),
        strip.background=element_rect(color="black"))+
  scale_fill_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                     values=cbPalette)+
  scale_shape_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                     values=c(21,22,24))+
  xlab("")+
  ylab(expression(paste("Net biodiversity effect (g ", m^-2, ")")))+
  geom_text(aes(x=-1, y=600, label="a"), 
            size=5, fontface="bold", inherit.aes = FALSE)

model<-data.frame(
  newx=seq(min(L2013$P),max(L2013$P),by=0.01),
  prd=modelCE$coef$`2013`$`coef(SMA)`[2]*seq(min(L2013$P),max(L2013$P),by=0.01)+modelCE$coef$`2013`$`coef(SMA)`[1])

p2<-ggplot(L2013, aes(x=P, y=CE, fill=Arrival, shape=Arrival))+
  geom_segment(aes(x=CIlow, y=CE, xend=CIhigh, yend=CE), color="gray50")+
  geom_point(color="black", size=2)+
  geom_line(aes(x=newx, y=prd), data=model, color="black", size=1, inherit.aes = FALSE)+
  theme_bw()+
  scale_y_continuous(breaks=seq(-200,500, by=100), limits=c(-200,500))+
  geom_vline(xintercept = 0, linetype=3, color="black")+
  geom_hline(yintercept = 0, linetype=3, color="black")+
  theme(legend.position="none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size=cex.axis,color="black"),
        axis.text.y=element_text(size=cex.axis,color="black"),
        axis.title.y=element_text(size=cex.title, margin=margin(t=0,r=7,b=0,l=0)),
        axis.title.x=element_text(size=cex.title, margin=margin(t=7,r=0,b=0,l=0)),
        strip.text=element_text(size=13),
        panel.border = element_rect(color="black"),
        strip.background=element_rect(color="black"))+
  scale_fill_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                    values=cbPalette)+
  scale_shape_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                     values=c(21,22,24))+
  xlab("")+
  ylab(expression(paste("Complementarity effect (g ", m^-2, ")")))+
  geom_text(aes(x=-1, y=500, label="c"), 
            size=5, fontface="bold", inherit.aes = FALSE)

p3<-ggplot(L2013, aes(x=P, y=DE, fill=Arrival, shape=Arrival))+
  geom_segment(aes(x=CIlow, y=DE, xend=CIhigh, yend=DE), color="gray50")+
  geom_point(color="black", size=2)+
  theme_bw()+
  scale_y_continuous(breaks=seq(0,300, by=50), limits=c(0,300))+
  geom_vline(xintercept = 0, linetype=3, color="black")+
  theme(legend.position="none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size=cex.axis,color="black"),
        axis.text.y=element_text(size=cex.axis,color="black"),
        axis.title.y=element_text(size=cex.title, margin=margin(t=0,r=7,b=0,l=0)),
        axis.title.x=element_text(size=cex.title, margin=margin(t=7,r=0,b=0,l=0)),
        strip.text=element_text(size=13),
        panel.border = element_rect(color="black"),
        strip.background=element_rect(color="black"))+
  scale_fill_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                    values=cbPalette)+
  scale_shape_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                     values=c(21,22,24))+
  xlab(expression(paste("Priority effect index (", italic(P), ")", sep="")))+
  ylab(expression(paste("Dominance effect (g ", m^-2, ")")))+
  geom_text(aes(x=0.72, y=300, label="e"), 
            size=5, fontface="bold", inherit.aes = FALSE)

model<-data.frame(
  newx=seq(min(L2014$P),max(L2014$P),by=0.01),
  prd=modelNBE$coef$`2014`$`coef(SMA)`[2]*seq(min(L2014$P),max(L2014$P),by=0.01)+modelNBE$coef$`2014`$`coef(SMA)`[1])

p4<-ggplot(L2014, aes(x=P, y=NBE, fill=Arrival, shape=Arrival))+
  geom_segment(aes(x=CIlow, y=NBE, xend=CIhigh, yend=NBE), color="gray50")+
  geom_point(color="black", size=2)+
  geom_line(aes(x=newx, y=prd), data=model, color="black", size=1, inherit.aes = FALSE)+
  theme_bw()+
  facet_wrap(~Year)+
  scale_y_continuous(breaks=seq(-100,300, by=100), limits=c(-100,300))+
  geom_vline(xintercept = 0, linetype=3, color="black")+
  geom_hline(yintercept = 0, linetype=3, color="black")+
  theme(legend.position="none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size=cex.axis,color="black"),
        axis.text.y=element_text(size=cex.axis,color="black"),
        axis.title.y=element_text(size=cex.title, margin=margin(t=0,r=7,b=0,l=0)),
        axis.title.x=element_text(size=cex.title, margin=margin(t=7,r=0,b=0,l=0)),
        strip.text=element_text(size=12),
        panel.border = element_rect(color="black"),
        strip.background=element_rect(color="black"))+
  scale_fill_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                    values=cbPalette)+
  scale_shape_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                     values=c(21,22,24))+
  xlab("")+
  ylab("")+
  geom_text(aes(x=-1, y=290, label="b"), 
            size=5, fontface="bold", inherit.aes = FALSE)

model<-data.frame(
  newx=seq(min(L2014$P),max(L2014$P),by=0.01),
  prd=modelCE$coef$`2014`$`coef(SMA)`[2]*seq(min(L2014$P),max(L2014$P),by=0.01)+modelCE$coef$`2014`$`coef(SMA)`[1])

p5<-ggplot(L2014, aes(x=P, y=CE, fill=Arrival, shape=Arrival))+
  geom_segment(aes(x=CIlow, y=CE, xend=CIhigh, yend=CE), color="gray50")+
  geom_point(color="black", size=2)+
  geom_line(aes(x=newx, y=prd), data=model, color="black", size=1, inherit.aes = FALSE)+
  theme_bw()+
  scale_y_continuous(breaks=seq(-100,300, by=100), limits=c(-100,300))+
  geom_vline(xintercept = 0, linetype=3, color="black")+
  geom_hline(yintercept = 0, linetype=3, color="black")+
  theme(legend.position="none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size=cex.axis,color="black"),
        axis.text.y=element_text(size=cex.axis,color="black"),
        axis.title.y=element_text(size=cex.title, margin=margin(t=0,r=7,b=0,l=0)),
        axis.title.x=element_text(size=cex.title, margin=margin(t=7,r=0,b=0,l=0)),
        strip.text=element_text(size=13),
        panel.border = element_rect(color="black"),
        strip.background=element_rect(color="black"))+
  scale_fill_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                    values=cbPalette)+
  scale_shape_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                     values=c(21,22,24))+
  xlab("")+
  ylab("")+
  geom_text(aes(x=-1, y=290, label="d"), 
            size=5, fontface="bold", inherit.aes = FALSE)

p6<-ggplot(L2014, aes(x=P, y=DE, fill=Arrival, shape=Arrival))+
  geom_segment(aes(x=CIlow, y=DE, xend=CIhigh, yend=DE), color="gray50")+
  geom_point(color="black", size=2)+
  theme_bw()+
  scale_y_continuous(breaks=seq(-100,150, by=50), limits=c(-100,150))+
  geom_vline(xintercept = 0, linetype=3, color="black")+
  geom_hline(yintercept = 0, linetype=3, color="black")+
  theme(legend.position="none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size=cex.axis,color="black"),
        axis.text.y=element_text(size=cex.axis,color="black"),
        axis.title.y=element_text(size=cex.title, margin=margin(t=0,r=7,b=0,l=0)),
        axis.title.x=element_text(size=cex.title, margin=margin(t=7,r=0,b=0,l=0)),
        strip.text=element_text(size=13),
        panel.border = element_rect(color="black"),
        strip.background=element_rect(color="black"))+
  scale_fill_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                    values=cbPalette)+
  scale_shape_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                     values=c(21,22,24))+
  xlab(expression(paste("Priority effect index (", italic(P), ")", sep="")))+
  ylab("")+
  geom_text(aes(x=-1, y=145, label="f"), 
            size=5, fontface="bold", inherit.aes = FALSE)

p7<-ggarrange(p1,p4,p2,p5,p3,p6, nrow=3, ncol=2, heights=c(1.1,1,1))

ggsave("Figure3.tiff", plot=p7, dpi=600, compression="lzw", 
       width=16, height=20, units="cm", pointsize=10)

###########################

plot(P2013, P2014, ylab="Priority effect index 2014", xlab="Priority effect index 2013",
     pch=c(rep(21,4), rep(22,4), rep(24,4)), las=1)
abline(b=1, a=0)

#############################
#Plot supplementary figure 4
#############################

cbPalette <- c("#0072B2", "#F0E442", "#009E73")

cex.axis=10
cex.title=10

model<-data.frame(
  newx=seq(min(L2013$P),max(L2013$P),by=0.01),
  prd=modelTICE$coef$`2013`$`coef(SMA)`[2]*seq(min(L2013$P),max(L2013$P),by=0.01)+modelTICE$coef$`2013`$`coef(SMA)`[1])

p1<-ggplot(L2013, aes(x=P, y=TICE, fill=Arrival, shape=Arrival))+
  geom_segment(aes(x=CIlow, y=TICE, xend=CIhigh, yend=TICE), color="gray50")+
  geom_point(color="black", size=2)+
  geom_line(aes(x=newx, y=prd), data=model, color="black", size=1, inherit.aes = FALSE)+
  theme_bw()+
  facet_wrap(~Year)+
  scale_y_continuous(breaks=seq(-200,400, by=100), limits=c(-200,400))+
  geom_vline(xintercept = 0, linetype=3, color="black")+
  geom_hline(yintercept = 0, linetype=3, color="black")+
  theme(legend.title=element_blank(),
        legend.key.size=unit(0.8, units="lines"),
        legend.position=c(0.85,0.18),
        legend.background = element_rect(fill=NA),
        panel.grid=element_blank(),
        axis.text.x=element_text(size=cex.axis,color="black"),
        axis.text.y=element_text(size=cex.axis,color="black"),
        axis.title.y=element_text(size=cex.title, margin=margin(t=0,r=7,b=0,l=0)),
        axis.title.x=element_text(size=cex.title, margin=margin(t=7,r=0,b=0,l=0)),
        strip.text=element_text(size=12),
        panel.border = element_rect(color="black"),
        strip.background=element_rect(color="black"))+
  scale_fill_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                    values=cbPalette)+
  scale_shape_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                     values=c(21,22,24))+
  xlab("")+
  ylab(expression(paste("TICE (g ", m^-2, ")")))+
  geom_text(aes(x=-1, y=400, label="a"), 
            size=5, fontface="bold", inherit.aes = FALSE)

model<-data.frame(
  newx=seq(min(L2013$P),max(L2013$P),by=0.01),
  prd=modelTDCE$coef$`2013`$`coef(SMA)`[2]*seq(min(L2013$P),max(L2013$P),by=0.01)+modelTDCE$coef$`2013`$`coef(SMA)`[1])

p2<-ggplot(L2013, aes(x=P, y=TDCE, fill=Arrival, shape=Arrival))+
  geom_segment(aes(x=CIlow, y=TDCE, xend=CIhigh, yend=TDCE), color="gray50")+
  geom_point(color="black", size=2)+
  geom_line(aes(x=newx, y=prd), data=model, color="black", size=1, inherit.aes = FALSE)+
  theme_bw()+
  scale_y_continuous(breaks=seq(-100,150, by=50), limits=c(-100,150))+
  geom_vline(xintercept = 0, linetype=3, color="black")+
  geom_hline(yintercept = 0, linetype=3, color="black")+
  theme(legend.position="none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size=cex.axis,color="black"),
        axis.text.y=element_text(size=cex.axis,color="black"),
        axis.title.y=element_text(size=cex.title, margin=margin(t=0,r=7,b=0,l=0)),
        axis.title.x=element_text(size=cex.title, margin=margin(t=7,r=0,b=0,l=0)),
        strip.text=element_text(size=13),
        panel.border = element_rect(color="black"),
        strip.background=element_rect(color="black"))+
  scale_fill_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                    values=cbPalette)+
  scale_shape_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                     values=c(21,22,24))+
  xlab(expression(paste("Priority effect index (", italic(P), ")", sep="")))+
  ylab(expression(paste("TDCE (g ", m^-2, ")")))+
  geom_text(aes(x=-1, y=150, label="c"), 
            size=5, fontface="bold", inherit.aes = FALSE)

model<-data.frame(
  newx=seq(min(L2014$P),max(L2014$P),by=0.01),
  prd=modelTICE$coef$`2014`$`coef(SMA)`[2]*seq(min(L2014$P),max(L2014$P),by=0.01)+modelTICE$coef$`2014`$`coef(SMA)`[1])

p3<-ggplot(L2014, aes(x=P, y=TICE, fill=Arrival, shape=Arrival))+
  geom_segment(aes(x=CIlow, y=TICE, xend=CIhigh, yend=TICE), color="gray50")+
  geom_point(color="black", size=2)+
  geom_line(aes(x=newx, y=prd), data=model, color="black", size=1, inherit.aes = FALSE)+
  theme_bw()+
  facet_wrap(~Year)+
  scale_y_continuous(breaks=seq(-100,300, by=100), limits=c(-100,300))+
  geom_vline(xintercept = 0, linetype=3, color="black")+
  geom_hline(yintercept = 0, linetype=3, color="black")+
  theme(legend.position="none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size=cex.axis,color="black"),
        axis.text.y=element_text(size=cex.axis,color="black"),
        axis.title.y=element_text(size=cex.title, margin=margin(t=0,r=7,b=0,l=0)),
        axis.title.x=element_text(size=cex.title, margin=margin(t=7,r=0,b=0,l=0)),
        strip.text=element_text(size=12),
        panel.border = element_rect(color="black"),
        strip.background=element_rect(color="black"))+
  scale_fill_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                    values=cbPalette)+
  scale_shape_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                     values=c(21,22,24))+
  xlab("")+
  ylab("")+
  geom_text(aes(x=-1, y=290, label="b"), 
            size=5, fontface="bold", inherit.aes = FALSE)

model<-data.frame(
  newx=seq(min(L2014$P),max(L2014$P),by=0.01),
  prd=modelTDCE$coef$`2014`$`coef(SMA)`[2]*seq(min(L2014$P),max(L2014$P),by=0.01)+modelTDCE$coef$`2014`$`coef(SMA)`[1])

p4<-ggplot(L2014, aes(x=P, y=TDCE, fill=Arrival, shape=Arrival))+
  geom_segment(aes(x=CIlow, y=TDCE, xend=CIhigh, yend=TDCE), color="gray50")+
  geom_point(color="black", size=2)+
  geom_line(aes(x=newx, y=prd), data=model, color="black", size=1, inherit.aes = FALSE)+
  theme_bw()+
  scale_y_continuous(breaks=seq(-40,40, by=20), limits=c(-40,40))+
  geom_vline(xintercept = 0, linetype=3, color="black")+
  geom_hline(yintercept = 0, linetype=3, color="black")+
  theme(legend.position="none",
        panel.grid=element_blank(),
        axis.text.x=element_text(size=cex.axis,color="black"),
        axis.text.y=element_text(size=cex.axis,color="black"),
        axis.title.y=element_text(size=cex.title, margin=margin(t=0,r=7,b=0,l=0)),
        axis.title.x=element_text(size=cex.title, margin=margin(t=7,r=0,b=0,l=0)),
        strip.text=element_text(size=13),
        panel.border = element_rect(color="black"),
        strip.background=element_rect(color="black"))+
  scale_fill_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                    values=cbPalette)+
  scale_shape_manual(name="Arrival", labels=c("F-first", "G-first", "L-first"), 
                     values=c(21,22,24))+
  xlab(expression(paste("Priority effect index (", italic(P), ")", sep="")))+
  ylab("")+
  geom_text(aes(x=-1, y=38, label="d"), 
            size=5, fontface="bold", inherit.aes = FALSE)

p5<-ggarrange(p1,p3,p2,p4, nrow=2, ncol=2, heights=c(1.1,1))

ggsave("FigureS4.tiff", plot=p5, dpi=600, compression="lzw", 
       width=16, height=13.5, units="cm", pointsize=8)
