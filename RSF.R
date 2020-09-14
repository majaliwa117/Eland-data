####################################################################
##################   Resource selection function ###################
# In this script I intend to answer one research question:
#i) Do eland preferentially select areas with high resources abundance or do they select areas free from natural predation across seasons? 

#clearing the R working environment
rm(list=ls())

######################Loading packages #############################
rm(list=ls())
options(warn=-1)
require(lme4)
require(MASS)
require(AICcmodavg)
require(performance)#Check for multicollinearity 
#check_collinearity()
#multicollinearity()

#Setting working directory
setwd("C:\\Users\\majal\\Documents\\MSc analysis")

#Listing files
list.files()

#Importing files to R
eland<- read.csv("RSF_data.csv", header=T)

attach(eland)#Attach eland

names(eland)#Gives the names of the column

#Checking correlation of explanatory variables
explvars<-cbind(D_forbs,Vegetation.cover,Vegetation.biomass,Heap.dung, 
                Escape.impediment,Rainfall,NDVI,dNDVI,dist_rangerpost,
                dist_confluence,dist_road,dist_drainage,dist_river,dist_parkbound,
                Elevation,slope,TWI,Curvature,Tourism_footprint,Human_proximity,Tree_density)
  
cor.explvars<-cor(explvars, use="complete.obs")# correlation matrix
cor.explvars<- data.frame(cor.explvars)#putting in a dataframe
write.csv(cor.explvars,"correlation.csv")#save it to CsV
#Elevation is left out because of strong correlation with other variables

############################ Standardizing the data #################
eland$D_forbs<-(D_forbs-mean(D_forbs))/sd(D_forbs)
eland$Vegetation.cover<-(Vegetation.cover-mean(Vegetation.cover))/sd(Vegetation.cover)
eland$Vegetation.biomass<-(Vegetation.biomass-mean(Vegetation.biomass))/sd(Vegetation.biomass)
eland$dist_rangerpost<-(dist_rangerpost-mean(dist_rangerpost))/sd(dist_rangerpost)
eland$dist_road<-(dist_road-mean(dist_road))/sd(dist_road)
eland$dist_river<-(dist_river-mean(dist_river))/sd(dist_river)
eland$dist_parkbound<-(dist_parkbound-mean(dist_parkbound))/sd(dist_parkbound)
eland$slope<-(slope-mean(slope))/sd(slope)
eland$Curvature<-(Curvature-mean(Curvature))/sd(Curvature)
eland$TWI<-(TWI-mean(TWI))/sd(TWI)
eland$Tourism_footprint<- (Tourism_footprint-mean(Tourism_footprint))/sd(Tourism_footprint)
eland$Heap.dung<-(Heap.dung-mean(Heap.dung))/sd(Heap.dung)
eland$Escape.impediment<-(Escape.impediment-mean(Escape.impediment))/sd(Escape.impediment)
eland$Human_proximity<-(Human_proximity-mean(Human_proximity))/sd(Human_proximity)
eland$Tree_density<-(Tree_density-mean(Tree_density,na.rm=T))/sd(Tree_density,na.rm=T)
eland$Rainfall<-(Rainfall-mean(Rainfall,na.rm=T))/sd(Rainfall,na.rm=T)
eland$NDVI<-(NDVI-mean(NDVI,na.rm=T))/sd(NDVI,na.rm=T)
eland$dNDVI<-(dNDVI-mean(dNDVI,na.rm=T))/sd(dNDVI,na.rm=T)
eland$dist_confluence<-(dist_confluence-mean(dist_confluence,na.rm=T))/sd(dist_confluence,na.rm=T)
eland$dist_drainage<-(dist_drainage-mean(dist_drainage,na.rm=T))/sd(dist_drainage,na.rm=T)
detach(eland)
attach(eland)
###################################################################
###################### Dividing the data seasonally################
elandD<-subset(eland,Season=="Dry")#Dry season
elandW<-subset(eland, Season=="Wet")#Wet season

detach(eland)#Detaching eland data
attach(elandD)#Attaching dry season

##################### Modelling dry season ######################
##Global model
Md1<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + dist_rangerpost + Tree_density+ 
             dist_road + Rainfall + NDVI +dNDVI +dist_parkbound + 
             Human_proximity + Tourism_footprint + Micro.habitat.type + 
             dist_confluence + dist_drainage +slope +
             Curvature+ slope +TWI+ dist_river+ Escape.impediment+
             I(dist_rangerpost^2) + I(dist_road^2) + I(dist_parkbound^2)+
             I(Human_proximity^2) +I(dist_confluence^2)+ I(dist_drainage^2) +
             I(dist_river^2) + (1|Animal_ID), family=binomial)

##M2(Remove distance to river, confluence)
Md2<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + Heap.dung + dist_rangerpost + Tree_density+ 
             dist_road + Rainfall + NDVI +dNDVI + dist_parkbound + Human_proximity +Tourism_footprint + 
            Micro.habitat.type + dist_drainage +slope+ Curvature+ I(NDVI^2)+ I(dist_parkbound^2)+
            I(dist_drainage^2)+I(dist_rangerpost^2) + I(Human_proximity^2)+ 
            (1|Animal_ID), family=binomial)

#M3(Remove distnace to drainage,confluence)
Md3<-glmer(Use~D_forbs + Vegetation.biomass + Heap.dung + Tree_density +
            NDVI +dNDVI + slope + Rainfall +Vegetation.cover + dist_river + 
            Escape.impediment+Curvature +dist_parkbound + dist_rangerpost + 
            dist_road + Tourism_footprint + Human_proximity+ Micro.habitat.type + 
            Vegetation.cover*dist_river + I(NDVI^2) + I(dist_river^2)+
            I(dist_parkbound^2)+ I(dist_rangerpost^2)+ I(dist_road^2)+
            I(Distance_humans^2)+ (1|Animal_ID), family=binomial)

#Remove distance to river,drainage,and confluence
Md4<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + slope +
             dist_rangerpost + dist_road + dist_parkbound + dist_river +
             Tourism_footprint + Micro.habitat.type + Human_proximity + Tree_density + 
             NDVI +dNDVI + Rainfall + Rainfall*(D_forbs + Vegetation.cover + 
             Vegetation.biomass+ dNDVI)+ Vegetation.biomass*Human_proximity +
             (1|Animal_ID), family=binomial)

##Remove distance to river,drainage,and fence
Md5<- glmer(Use~D_forbs + Vegetation.biomass + Heap.dung + Curvature+ Escape.impediment+
             dist_rangerpost + dist_road + dist_parkbound + Tourism_footprint + Micro.habitat.type +
             dist_confluence + slope + Human_proximity+ Tree_density + Rainfall +NDVI +dNDVI +
             Human_proximity*Vegetation.biomass + I(NDVI^2) + I(dist_confluence^2)+
             I(dist_parkbound^2)+ I(dist_rangerpost^2)+ I(dist_road^2)+
             I(Human_proximity^2)+(1|Animal_ID), family=binomial)

#Remove distance to confluence,TWI, and vegetation biomass
Md6<-glmer(Use~D_forbs +  Heap.dung + Tree_density + dist_river + Vegetation.cover + 
            dist_rangerpost + dist_road + dist_parkbound +Tourism_footprint + 
            slope + Rainfall + NDVI +dNDVI +dist_drainage +Micro.habitat.type +
            I(NDVI^2) + I(Proximity_river^2)+ I(dist_parkbound^2)+ I(dist_rangerpost^2)+
            I(dist_road^2)+ I(Human_proximity^2)+I(dist_drainage^2)+
            (1|Animal_ID),family = binomial)

best_modelD<-list(Md1,Md2,Md3,Md4,Md5,Md6)# selecting the best model

print(as.data.frame(aictab(best_modelD)),sort=T)# 
##Md4 has lowest AIC and highest weight

## Checking collinearity
multicollinearity(Md4)

sink("Dry.txt") #saving the results in text
print(summary(Md4))#Getting summary
sink()#saving it instead of printing
###################################################################
############################# Wet season Modelling#################
detach(elandD)
attach(elandW)
##Global model
Mw1<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + 
             Heap.dung + dist_rangerpost + Tree_density+ 
             dist_road + Rainfall + NDVI +dNDVI + dist_parkbound + 
             Human_proximity + Tourism_footprint + Micro.habitat.type + 
             dist_confluence + dist_drainage +slope + Curvature+ slope + 
             TWI + dist_river+ Escape.impediment+I(NDVI^2)+ I(dist_rangerpost^2) + 
             I(dist_road^2) + I(dist_confluence^2)+ I(dist_parkbound^2) +
             I(Human_proximity^2)+ I(dist_drainage^2)+
             I(dist_river^2)+(1|Animal_ID), family=binomial)

#M2(Remove distance drainage,confluence)
Mw2<-glmer(Use~D_forbs + Vegetation.biomass + Heap.dung + Tree_density +
             NDVI +dNDVI + slope + Rainfall +Vegetation.cover +dist_river + 
             Escape.impediment+ Curvature + dist_parkbound + dist_rangerpost + 
             dist_road + Tourism_footprint + Human_proximity + Micro.habitat.type + 
             Vegetation.cover*dist_river +I(NDVI^2)+ I(dist_rangerpost^2) + I(dist_road^2) + 
             I(dist_parkbound^2) +I(Human_proximity^2)+ I(dist_river^2)+
             (1|Animal_ID), family=binomial)

#Remove distance to river,drainage and confluence
Mw3<-glmer(Use~D_forbs + Vegetation.cover + Vegetation.biomass + Heap.dung + Curvature +
             Escape.impediment+ dist_rangerpost + dist_road + dist_parkbound +
             Tourism_footprint + Micro.habitat.type + slope + Human_proximity+
             Tree_density + Rainfall +NDVI +dNDVI + Rainfall*(D_forbs + Vegetation.cover + 
             Vegetation.biomass + NDVI + dNDVI)+ Vegetation.biomass*Human_proximity +I(NDVI^2)+
             I(dist_rangerpost^2) + I(dist_road^2) + I(dist_parkbound^2) +I(Human_proximity^2)+
            (1|Animal_ID), family=binomial)

##Remove distance to river and drainage
Mw4<- glmer(Use~D_forbs + Vegetation.biomass + Heap.dung + Curvature+ Escape.impediment+
              dist_rangerpost + dist_road + dist_parkbound +
              Tourism_footprint + Micro.habitat.type + dist_confluence + 
              slope + Human_proximity+ Tree_density + Rainfall +NDVI +dNDVI +
              Human_proximity*Vegetation.biomass +I(NDVI^2)+ I(dist_road^2)+ I(dist_rangerpost^2) +
              I(dist_road^2)+ I(dist_parkbound^2) +I(dist_confluence^2)+
              (1|Animal_ID), family=binomial)

#Remove proximity to confluence,TWI, and vegetation biomass
Mw5<-glmer(Use~D_forbs +Vegetation.biomass+  Heap.dung + Tree_density + dist_river + Vegetation.cover + 
             dist_rangerpost + dist_road + dist_parkbound + Tourism_footprint + slope +
             Rainfall + NDVI +dNDVI + dist_drainage + Micro.habitat.type +I(NDVI^2)+
             I(dist_river^2)+ I(dist_rangerpost^2) + I(dist_road^2)+
             I(dist_parkbound^2) +I(dist_drainage^2)+(1|Animal_ID),family = binomial)

#
Mw6<-glmer(Use~D_forbs +Vegetation.biomass + TWI + Curvature +
             Rainfall + NDVI + dNDVI+ Vegetation.cover+ dist_parkbound +
             Tourism_footprint + dist_rangerpost + dist_road + 
             Tree_density + Micro.habitat.type + I(NDVI^2)+ I(dist_parkbound^2)+ 
             (1|Animal_ID),family = binomial)

#Model selection
best_modelW<-list(Mw1,Mw2,Mw3,Mw4,Mw5,Mw6)# selecting the best model

#Printing the best model
print(as.data.frame(aictab(best_modelW)),sort=T)##Mw6 has lowest AIC 

## Checking collinearity
multicollinearity(Mw6)

sink("Wet.txt") #saving the results in text
print(summary(Mw6))#Getting summary
sink()#saving it instead of printing
###################################################################