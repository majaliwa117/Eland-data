##################################################################
##           Integrated step selection analysis ##################
##################################################################
# In this script I intend determine to answer two questions:
#(i)How do resource abundance, predation risk and anthropogenic disturbance influence how eland moves ? and 
#(ii)Does eland move faster or slower if they start in anthropogenic or predation risk area? 

#clearing the R working environment
rm(list=ls())

##Loading packages
require(survival)
require(AICcmodavg)
require(performance)#Checking multicollinearity

#Loading data
eland<-read.csv("C:\\Users\\majal\\Documents\\MSc analysis\\ISSF\\Issa_data.csv")

eland$tod_end_<-ifelse(eland$tod_end_=="day","Day","Night")#change the first letter to capital

attach(eland) #attaching data

#Checking correlation of explanatory variables
explvars<-cbind(log_sl,cos_ta,end_dist_confluence,end_dist_road,
                end_dist_river,end_Tourism_footprint,end_dist_rangerpost,
                end_dist_parkbound,end_Human_proximity,end_tree_density,
                end_rainfall,end_Elevation,end_TWI,end_Curvature,end_NDVI,end_dNDVI)

cor.explvars<-cor(explvars, use="complete.obs")# correlation matrix
cor.explvars<- data.frame(cor.explvars)#putting in a dataframe
#distance to river and distance to rangerpost as they are strongly
#related,they will not be included in the same model
#Drop elevation and proximity to confluence_2 because it save the same purpose as landscape curvature

## Standardizing the variables
#movement(start point)
eland$start_dist_road<-(start_dist_road-mean(start_dist_road))/sd(start_dist_road)
eland$start_dist_river<-(start_dist_river-mean(start_dist_river))/sd(start_dist_river)
eland$start_Tourism_footprint<-(start_Tourism_footprint-mean(start_Tourism_footprint))/sd(start_Tourism_footprint)
eland$start_dist_rangerpost<-(start_dist_rangerpost-mean(start_dist_rangerpost))/sd(start_dist_rangerpost)
eland$start_dist_parkbound<-(start_dist_parkbound-mean(start_dist_parkbound,na.rm=T))/sd(start_dist_parkbound,na.rm=T)
eland$start_dist_confluence<-(start_dist_confluence-mean(start_dist_confluence,na.rm=T))/sd(start_dist_confluence,na.rm=T)
eland$start_Human_proximity<-(start_Human_proximity-mean(start_Human_proximity,na.rm=T))/sd(start_Human_proximity,na.rm=T)
eland$start_tree_density<-(start_tree_density-mean(start_tree_density,na.rm=T))/sd(start_tree_density,na.rm=T)

#Resource selection(end point)
eland$end_TWI<-(end_TWI-mean(end_TWI,na.rm=T))/sd(end_TWI,na.rm=T)
eland$end_Curvature<-(end_Curvature-mean(end_Curvature))/sd(end_Curvature)
eland$end_rainfall<-(end_rainfall-mean(end_rainfall))/sd(end_rainfall)
eland$end_NDVI<-(end_NDVI-mean(end_NDVI,na.rm=T))/sd(end_NDVI,na.rm=T)
eland$end_dNDVI<-(end_dNDVI-mean(end_dNDVI,na.rm=T))/sd(end_dNDVI,na.rm=T)
eland$end_dist_road<-(end_dist_road-mean(end_dist_road))/sd(end_dist_road)
eland$end_dist_river<-(end_dist_river-mean(end_dist_river))/sd(end_dist_river)
eland$end_Tourism_footprint<-(end_Tourism_footprint-mean(end_Tourism_footprint,na.rm=T))/sd(end_Tourism_footprint,na.rm=T)
eland$end_dist_rangerpost<-(end_dist_rangerpost-mean(end_dist_rangerpost))/sd(end_dist_rangerpost)
eland$end_dist_parkbound<-(end_dist_parkbound-mean(end_dist_parkbound,na.rm=T))/sd(end_dist_parkbound,na.rm=T)
eland$end_dist_confluence<-(end_dist_confluence-mean(end_dist_confluence,na.rm=T))/sd(end_dist_confluence,na.rm=T)
eland$end_Human_proximity<-(end_Human_proximity-mean(end_Human_proximity,na.rm=T))/sd(end_Human_proximity,na.rm=T)
eland$end_tree_density<-(end_tree_density-mean(end_tree_density,na.rm=T))/sd(end_tree_density,na.rm=T)

detach(eland)
################fitting model using clogit function #################
attach(eland)
#Null model
M1<-clogit(case_~cos_ta + log_sl + log_sl:cos_ta + strata(step_id_))

#Global model
M2<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
            end_TWI+ end_rainfall+ end_NDVI + end_dNDVI + end_dist_road + end_dist_river +
            end_Tourism_footprint + end_dist_rangerpost + end_dist_parkbound + end_Human_proximity +
            end_tree_density + I(end_NDVI^2)+ log_sl*(start_Tourism_footprint + start_dist_river +
            start_dist_road + start_dist_rangerpost + start_dist_parkbound + start_Human_proximity +
            start_tree_density) + strata(step_id_))

#Remove distance to river(endpoint)
M3<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
             end_TWI+ end_rainfall + end_NDVI + end_dNDVI + end_dist_road + end_Tourism_footprint + 
             end_dist_rangerpost + end_dist_rangerpost + end_Human_proximity + end_tree_density +
             I(end_NDVI^2)+ log_sl*(start_Tourism_footprint +start_dist_road + start_dist_rangerpost +
             start_dist_parkbound + start_Human_proximity + start_tree_density) + strata(step_id_))

#Remove distance to rangerpost(endpoint)
M4<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
           end_TWI+ end_rainfall + end_NDVI + end_dNDVI + end_dist_road + end_dist_river +
             end_Tourism_footprint + end_dist_parkbound + end_Human_proximity + end_tree_density + 
             I(end_NDVI^2)+ log_sl*(start_Tourism_footprint + start_dist_river +start_dist_road +
             start_dist_parkbound + start_Human_proximity + start_tree_density) + strata(step_id_))

#Model selection
best_model<-list(M1,M2,M3,M4)# selecting the best model
h<-data.frame(print(as.data.frame(aictab(best_model)),sort=T))#The best model is M2

####################################################################
#             Likelihood ratio test 
####################################################################
# Identifying the most important variable in the model and dropping the 
#Insignificant ones.
M2_1<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
               end_TWI+ end_rainfall+ end_NDVI + end_dNDVI + end_dist_road + end_dist_river +
               end_Tourism_footprint + end_dist_rangerpost + end_dist_parkbound + end_Human_proximity +
               end_tree_density + I(end_NDVI^2)+ log_sl*(start_Tourism_footprint + start_dist_river +
               start_dist_road + start_dist_rangerpost + start_dist_parkbound + start_Human_proximity +
               start_tree_density) + strata(step_id_))#Maximal model
logLik(M2_1)#Likelihood ratio ('log Lik.' -18655.16 (df=31))

#Droppring NDVI self interaction
M2_2<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
         end_TWI+ end_rainfall+ end_NDVI + end_dNDVI + end_dist_road + end_dist_river +
         end_Tourism_footprint + end_dist_rangerpost + end_dist_parkbound + end_Human_proximity +
         end_tree_density + log_sl*(start_Tourism_footprint + start_dist_river +
         start_dist_road + start_dist_rangerpost + start_dist_parkbound + start_Human_proximity +
         start_tree_density) + strata(step_id_))
logLik(M2_2)#Likelihood ratio ('log Lik.' -18658.46 (df=30)))

#Make comparison of the model(M1,M2)
1-pchisq(2*(-18655.16--18658.46),1)
#p=0.01019788 this suggest retention of maximal model

#Dropping of dNDVI(end point)
M2_3<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
               end_TWI+ end_rainfall+ end_NDVI + end_dist_road + end_dist_river +
               end_Tourism_footprint + end_dist_rangerpost + end_dist_parkbound + end_Human_proximity +
               end_tree_density + I(end_NDVI^2)+ log_sl*(start_Tourism_footprint + start_dist_river +
               start_dist_road + start_dist_rangerpost + start_dist_parkbound + start_Human_proximity +
               start_tree_density) + strata(step_id_))
logLik(M2_3)#Likelihood ratio('log Lik.' -18655.34 (df=30)))

#Make comparison of the model(M1,M3)
1-pchisq(2*(-18655.16--18655.34),1)
#p=0.5485062 this suggest adoption of simpler model(M2_3)

#Dropping rainfall
M2_4<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
               end_TWI+ end_NDVI + end_dist_road + end_dist_river +
               end_Tourism_footprint + end_dist_rangerpost + end_dist_parkbound + end_Human_proximity +
               end_tree_density + I(end_NDVI^2)+ log_sl*(start_Tourism_footprint + start_dist_river +
               start_dist_road + start_dist_rangerpost + start_dist_parkbound + start_Human_proximity +
              start_tree_density) + strata(step_id_))
logLik(M2_4)#Likelihood ratio ('log Lik.' -18655.39 (df=29)))

#Make comparison of the model(M3,M4)
1-pchisq(2*(-18655.34--18655.39),1)
#p=0.7518296 this suggest adoption of simpler model(M2_4)

#Dropping TWI (end point)
M2_5<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
               end_NDVI + end_dist_road + end_dist_river + end_Tourism_footprint + 
               end_dist_rangerpost + end_dist_parkbound + end_Human_proximity +
               end_tree_density + I(end_NDVI^2)+ log_sl*(start_Tourism_footprint + 
               start_dist_river +start_dist_road + start_dist_rangerpost + start_dist_parkbound +
               start_Human_proximity + start_tree_density) + strata(step_id_))
logLik(M2_5)#Likelihood ratio ('log Lik.' -18661.48 (df=28)))

#Make comparison of the model(M4,M5)
1-pchisq(2*(-18655.39--18661.48),1)
#p=0.0004830464 this suggest retention of model M2_4

#Dropping curvature (end point)
M2_6<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
               end_TWI+ end_NDVI + end_dist_road + end_dist_river +
               end_Tourism_footprint + end_dist_rangerpost + end_dist_parkbound + end_Human_proximity +
               end_tree_density + I(end_NDVI^2)+ log_sl*(start_Tourism_footprint + start_dist_river +
               start_dist_road + start_dist_rangerpost + start_dist_parkbound + start_Human_proximity +
               start_tree_density) + strata(step_id_))
logLik(M2_6)#Likelihood ratio ('log Lik.'-18658.26 (df=28)))

#Make comparison of the model(M4,M6)
1-pchisq(2*(-18655.39--18658.26),1)
#p=0.01658279 this suggest retention of model M2_4

#Dropping distance to river(end point)
M2_7<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
               end_TWI+ end_NDVI + end_dist_road + end_Tourism_footprint + end_dist_rangerpost + 
               end_dist_parkbound + end_Human_proximity + end_tree_density + I(end_NDVI^2)+ 
               log_sl*(start_Tourism_footprint + start_dist_river + start_dist_road + 
               start_dist_rangerpost + start_dist_parkbound + start_Human_proximity +
               start_tree_density) + strata(step_id_))
logLik(M2_7)#Likelihood ratio ('log Lik.'-18656.24 (df=28)))

#Make comparison of the model(M4,M7)
1-pchisq(2*(-18655.39--18656.24),1)
#p=0.192288 this suggest adoption of model M2_7

#Dropping Tourism footprint(end point)
M2_8<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
               end_TWI+ end_NDVI + end_dist_road + end_dist_rangerpost + 
               end_dist_parkbound + end_Human_proximity + end_tree_density + I(end_NDVI^2)+ 
               log_sl*(start_Tourism_footprint + start_dist_river + start_dist_road + 
               start_dist_rangerpost + start_dist_parkbound + start_Human_proximity +
               start_tree_density) + strata(step_id_))
logLik(M2_8)#Likelihood ratio ('log Lik.'-18656.76(df=27)))

#Make comparison of the model(M7,M8)
1-pchisq(2*(-18656.24--18656.76),1)
#p=0.3078215 this suggest adoption of model M2_8

#Dropping Human footprint(end point)
M2_9<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
               end_TWI+ end_NDVI + end_dist_road + end_dist_rangerpost + 
               end_dist_parkbound + end_tree_density + I(end_NDVI^2)+ 
               log_sl*(start_Tourism_footprint + start_dist_river + start_dist_road + 
               start_dist_rangerpost + start_dist_parkbound + start_Human_proximity +
               start_tree_density) + strata(step_id_))
logLik(M2_9)#Likelihood ratio('log Lik.'-18658.12(df=26)))

#Make comparison of the model(M8,M9)
1-pchisq(2*(-18656.76--18658.12),1)
#p=0.09909802 this suggest adoption of model M2_9

#Dropping Tourism footprint(start point)
M2_10<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
                end_TWI+ end_NDVI + end_dist_road + end_dist_rangerpost + 
                end_dist_parkbound + end_tree_density + I(end_NDVI^2)+ 
                log_sl*(start_dist_river + start_dist_road + start_dist_rangerpost + 
                start_dist_parkbound + start_Human_proximity + start_tree_density) + strata(step_id_))
logLik(M2_10)#Likelihood ratio('log Lik.'-18658.92(df=24)))

#Make comparison of the model(M9,M10)
1-pchisq(2*(-18658.12--18658.92),2)
#p= 0.449329 this suggest adoption of model M2_10

#Dropping distance to river(start point)
M2_11<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
              end_TWI+ end_NDVI + end_dist_road + end_dist_rangerpost + 
              end_dist_parkbound + end_tree_density + I(end_NDVI^2)+ 
              log_sl*(start_dist_road + start_dist_rangerpost + start_dist_parkbound +
              start_Human_proximity + start_tree_density) + strata(step_id_))
logLik(M2_11)#Likelihood ratio('log Lik.'-18665.74(df=22)))

#Make comparison of the model(M9,M10)
1-pchisq(2*(-18658.92--18665.74),2)
#p= 0.001091721 this suggest retention of model M2_10

#Dropping distance to park boundary(start point)
M2_12<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
                end_TWI+ end_NDVI + end_dist_road + end_dist_rangerpost + 
                end_dist_parkbound + end_tree_density + I(end_NDVI^2)+ 
                log_sl*(start_dist_river + start_dist_road + 
                start_dist_rangerpost + start_Human_proximity +
                start_tree_density) + strata(step_id_))
logLik(M2_12)#Likelihood ratio('log Lik.'-18675.04(df=22)))

#Make comparison of the model(M10,M12)
1-pchisq(2*(-18658.92--18675.04),2)
#p= 9.980975e-08 this suggest retention of  model M2_10

#Dropping human footprint(end point)
M2_13<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
                end_TWI+ end_NDVI + end_dist_road + end_dist_rangerpost + 
                end_dist_parkbound + end_tree_density + I(end_NDVI^2)+ 
                log_sl*(start_dist_river + start_dist_road + 
                start_dist_rangerpost + start_dist_parkbound + 
                start_tree_density) + strata(step_id_))
logLik(M2_13)#Likelihood ratio('log Lik.'-18659.18(df=22)))

#Make comparison of the model(M10,M13)
1-pchisq(2*(-18658.92--18659.18),2)
#p= 0.7710516 this suggest adoption of model M2_13

#Dropping distance to road(start point)
M2_14<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
                end_TWI+ end_NDVI + end_dist_road + end_dist_rangerpost + 
                end_dist_parkbound + end_tree_density + I(end_NDVI^2)+ 
                log_sl*(start_dist_river + start_dist_rangerpost + start_dist_parkbound + 
                start_tree_density) + strata(step_id_))
logLik(M2_14)#Likelihood ratio('log Lik.'-18668.62 (df=20)))

#Make comparison of the model(M13,M14)
1-pchisq(2*(-18659.18--18668.62 ),2)
#p= 7.948041e-05 this suggest retention of model M2_13

#Dropping distance to rangerpost(start point)
M2_15<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
                end_TWI+ end_NDVI + end_dist_road + end_dist_rangerpost + 
                end_dist_parkbound + end_tree_density + I(end_NDVI^2)+ 
                log_sl*(start_dist_river + start_dist_road + start_dist_parkbound + 
                start_tree_density) + strata(step_id_))
logLik(M2_15)#Likelihood ratio('log Lik.'-18661.88(df=20)))

#Make comparison of the model(M13,M15)
1-pchisq(2*(-18659.18--18661.88),2)
#p=0.06720551 this suggest adoption of model M2_15

#Summary of the best model
best_M<-clogit(case_~cos_ta + log_sl + cos_ta:log_sl + log_sl:tod_end_ + end_Curvature + 
                end_TWI+ end_NDVI + end_tree_density + end_dist_river+ end_dist_road + 
                end_dist_rangerpost + end_dist_parkbound + I(end_NDVI^2)+ 
                log_sl*(start_tree_density + start_dist_river + start_dist_road +
                start_dist_rangerpost + start_dist_parkbound) + strata(step_id_))

sink("bestmodel.txt") #saving the results in text
print(summary(best_M))#Getting summary
sink()#saving it instead of printing
####################################################################