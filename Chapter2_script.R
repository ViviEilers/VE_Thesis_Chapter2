#setwd("C:/Users/vivi_/Dropbox/PhD Vivi/GLM_RK") # Dell address
setwd("C:/Users/Vivianne Eilers/Dropbox/PhD Vivi\\GLM_RK") # Lenovo and HP address
#setwd("C:\\Users\\r01ve17\\Dropbox\\PhD Vivi\\GLM_RK") # Uni address

# Reading the data ####

#'*Model including all species in the dataset*
#'*Model gam.NB.random5b used for my first data chapter*

library(data.table)
library(mgcv)
library(RColorBrewer)
library(plyr)

RK7Rp.final<- fread("Data_RK_noZero.csv", data.table= FALSE, stringsAsFactors = TRUE) 

system.time(gam.NB.random5b<- bam(Quantity ~ -1 + LogBM + Car_velocity +
                                    Car_velocity:LogBM +
                                    Gap_Size + offset(LogSection) +
                                    s(SpBio, bs= "re") +
                                    s(S.Gs4plus, bs= "re"),
                                    data= RK7Rp.final, family= nb, 
                                    discrete= F, nthreads= 4)) 

save(gam.NB.random5b, file= "BAM5b_outputs.RData")
load("BAM5b_outputs.RData")
summary(gam.NB.random5b)

#'* summary(gam.NB.random5b) # from BAM5b_outputs.RData*
# 
# Family: Negative Binomial(1.065) 
# Link function: log 
# 
# Formula:
#   Quantity ~ -1 + LogBM + Car_velocity + Car_velocity:LogBM + Gap_Size + 
#   offset(LogSection) + s(SpBio, bs = "re") + s(S.Gs4plus, 
#                                                bs = "re")
# 
# Parametric coefficients:
#                        Estimate Std. Error t value Pr(>|t|)    
#   LogBM              -1.475716   0.083638 -17.644  < 2e-16 ***
#   Car_velocity       -0.337860   0.013908 -24.292  < 2e-16 ***
#   Gap_Size<4          3.208746   1.397095   2.297   0.0216 *  
#   Gap_Size>4          4.095056   0.604126   6.778 1.22e-11 ***
#   LogBM:Car_velocity  0.039242   0.001865  21.038  < 2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# Approximate significance of smooth terms:
#                  edf Ref.df     F p-value    
#   s(SpBio)     563.0    653 305.7  <2e-16 ***
#   s(S.Gs4plus) 144.3    454 437.2  <2e-16 ***
#   ---
#   Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# R-sq.(adj) =   0.18   Deviance explained = -78.1%
# fREML = 1.3508e+06  Scale est. = 1         n = 942509

#Model validation #####
par(mfrow=c(1,2),mar= c(4.1, 4.1, 2.1, 2.1), las=1)
plot(gam.NB.random5b)

par(mfrow=c(2,2))
gam.check(gam.NB.random5b)

#Extracting coefficients for all species ####
GNBR.coef5b<- coefficients(gam.NB.random5b)

par(mfrow=c(1,1))
hist(GNBR.coef5b[grep("SpBio", names(GNBR.coef5b))], nclass= 25)

SpBio.coefs5b<- GNBR.coef5b[grep("SpBio", names(GNBR.coef5b))] # It extracts the coefficient for all sp and biome

#Creating a data frame for the coefficient values and names
SpBio.df5b<- data.frame(coef= SpBio.coefs5b, code= names(SpBio.coefs5b))

#Adding variables to the data frame
SpBio5b<-unique(RK7Rp.final[, c("Fauna_group","Species", "Biome","adult_BM_g", "LogBM","IUCN","BRL")])
SpBio5b<- SpBio5b[order(SpBio5b$Biome),]
SpBio5b<- SpBio5b[order(SpBio5b$Species),]
SpBio.df5b$Sp<-SpBio5b$Species
SpBio.df5b$Bio<-SpBio5b$Biome
SpBio.df5b$LogBM<-SpBio5b$LogBM
SpBio.df5b$IUCN<-SpBio5b$IUCN
SpBio.df5b$BRL<-SpBio5b$BRL
SpBio.df5b$Fauna_group<-SpBio5b$Fauna_group
SpBio.df5b$adult_BM_g<-SpBio5b$adult_BM_g

#Calculating estimates for 100 km of the road per year for daily surveys at 40 km/h ####
SpBio.df5b$Estimate40<- exp(SpBio.df5b$coef +
                              40 * GNBR.coef5b["Car_velocity"] + # assuming 40 km/h survey
                              SpBio.df5b$LogBM * GNBR.coef5b["LogBM"] +
                              GNBR.coef5b["Gap_Size<4"] +
                              SpBio.df5b$LogBM * 40 * GNBR.coef5b["LogBM:Car_velocity"]) * 100 * 365 # unit is per 100 km per year

#Calculating the se.fit for the estimates at 40 km/h ####
SpBio.df5b$SpBio<-paste(SpBio.df5b$Sp, SpBio.df5b$Bio, sep="_")

Rk.40.5b<-data.frame(Car_velocity=40,
                     LogBM=SpBio.df5b$LogBM,
                     SpBio=SpBio.df5b$SpBio,
                     Gap_Size="<4",
                     LogSection=log(100*365),
                     S.Gs4plus="dummyLevel")


SE40_GNBR5b<-predict(gam.NB.random5b, newdata=Rk.40.5b, se.fit=T, type="response")
#str(SE40_GNBR5b)

SpBio.df5b$SE40<-SE40_GNBR5b$se.fit
SpBio.df5b$minSE40<-SpBio.df5b$Estimate40-SpBio.df5b$SE40
SpBio.df5b$maxSE40<-SpBio.df5b$Estimate40+SpBio.df5b$SE40 

#Calculating the relative standard error####
SpBio.df5b$PSE<-((SpBio.df5b$SE40 *100) / SpBio.df5b$Estimate40)

# Calculating sample size ####
SS.5b<-vector(mode = "logical", length = length(levels(RK7Rp.final$SpBio)))
for (x in 1: length(levels(RK7Rp.final$SpBio))){
  SS.5b[x]<-sum(subset(RK7Rp.final$Quantity,RK7Rp.final$SpBio==levels(RK7Rp.final$SpBio)[x]))
}

SS2.5b<-matrix(data = NA, nrow = length(levels(RK7Rp.final$SpBio)), ncol = 2)
SS2.5b[,1]<-levels(RK7Rp.final$SpBio)
SS2.5b[,2]<-SS.5b

SS2.5b<-as.data.frame(SS2.5b)
colnames(SS2.5b)<-c("SpBio","TotalCarcasses")

NewData.5b<-merge(RK7Rp.final, SS2.5b, by = "SpBio", all.x=TRUE, all.y=FALSE, no.dups=TRUE)
#summary(NewData.5b)
#head(NewData.5b)

TC<-unique(NewData.5b[, c("TotalCarcasses", "SpBio")])
SpBio.df5b$TotalCarcasses<-as.numeric(TC$TotalCarcasses)
SpBio.df5b$Sp.20<-sub("_", " ", SpBio.df5b$Sp)


##Plotting sample size against SE#####
par(mfrow=c(1,1),mar= c(4.1, 4.1, 3.1, 3.1), las=1)
plot(SpBio.df5b$TotalCarcasses, log(SpBio.df5b$PSE), xlab= "Sample size", ylab= "Standard error")
#Reducing the magnitude of differences
plot(log(SpBio.df5b$TotalCarcasses), log(SpBio.df5b$PSE), xlab= "Sample size (log)", ylab= "Standard error (log)") # Percentage of SE increases with sample size


# Exploring the data

# Calculating the total of carcasses per species and biome
#RK7Rp.final$TotalCarcasses<-NewData.5b$TotalCarcasses

#How many carcass species there are per biome? ####
summary(SpBio.df5b$Bio)

#In how many biomes each species is present? ####
summary(SpBio.df5b$Sp)

#How many species of big taxonomic group? #####
Sp.group5b<-unique(RK7Rp.final[, c("Species", "Fauna_group")])
summary(Sp.group5b)
#                 Species    Fauna_group
# Accipiter_striatus :  1   A:299      
# Agelaioides_badius :  1   H: 51      
# Agelasticus_thilius:  1   M:108      
# Akodon_cursor      :  1              
# Alouatta_caraya    :  1              
# Alouatta_guariba   :  1              
# (Other)            :452   

#How many carcasses of big taxonomic group? #####
C.group5b<-tapply(RK7Rp.final$Quantity, RK7Rp.final$Fauna_group, sum, na.rm=T)
C.group5b
#   A    H    M 
#5211 1589 6086

#How many species are in the IUCN category vulnerable?
summary(SpBio.df5b$Sp[SpBio.df5b$IUCN=="VU"])

#How many species are in the Brazilian red list category vulnerable?
summary(SpBio.df5b$Sp[SpBio.df5b$BRL=="VU"])

#How many species are in the Brazilian red list category endangered?
summary(SpBio.df5b$Sp[SpBio.df5b$BRL=="EN"])

#How many species are in the IUCN category endangered?
summary(SpBio.df5b$Sp[SpBio.df5b$IUCN=="EN"])

#Ranking the top 20 species recorded during the surveys####
RK.matrix.raw<- tapply(RK7Rp.final$Quantity, list(RK7Rp.final$Biome, RK7Rp.final$Species), sum)
RK.matrix.colmax.raw<- apply(RK.matrix.raw, MARGIN=2, FUN=max, na.rm= T)
RK.matrix.raw<- RK.matrix.raw[, order(RK.matrix.colmax.raw)]

barplot(RK.matrix.raw[, 439:458], beside=T, horiz=T, legend=rownames(RK.matrix.raw), 
        args.legend=list(x=1200, y=70, cex=0.8, text.width=120), col= brewer.pal(n=6, name="Set1"),
        xlab="Carcasses registered", main="Most abundant species")


# Ranking top 20 species estimated from model 5b by coefficient at 40 km/h ####
SpBio.matrix5b<- tapply(SpBio.df5b$Estimate40, list(SpBio.df5b$Bio, SpBio.df5b$Sp.20), mean)
SpBio.matrix.colmax5b<- apply(SpBio.matrix5b, MARGIN=2, FUN=max, na.rm= T)
SpBio.matrix5b<- SpBio.matrix5b[, order(SpBio.matrix.colmax5b)]
SpBio20.matrix5b<-SpBio.matrix5b[,439:458]

# Creating a matrix for SE40 for the top 20 species estimated at 40 km/h####
SpBio.matrix.SE5b<- tapply(SpBio.df5b$SE40, list(SpBio.df5b$Bio, SpBio.df5b$Sp.20), mean)
SpBio.matrix.SE5b<- SpBio.matrix.SE5b[, order(SpBio.matrix.colmax5b)]
SpBio.matrix.SE5b.20<-SpBio.matrix.SE5b[,439:458]

#Plotting the rank of the top 20 estimated 
par(mfrow=c(1,1), mar= c(4.1, 12, 3.1, 4.1), las=1)
Bp.Sp.5b<-barplot(SpBio20.matrix5b, beside=T, horiz=T, 
                  xlim=c(-12, 650), legend=rownames(SpBio20.matrix5b), 
                  args.legend=list(x="bottomright",cex=0.8,   
                                   text.width=80, bty="n"), col= brewer.pal(n=6, name="Set1"), font=3,
                  xlab="Estimated detection per 100 km/year", 
)
title(ylab="Species", line=11, cex.lab=1)

segments(SpBio20.matrix5b-SpBio.matrix.SE5b.20, Bp.Sp.5b,  
         SpBio20.matrix5b+SpBio.matrix.SE5b.20, Bp.Sp.5b,
         lwd=1.5)

arrows(SpBio20.matrix5b-SpBio.matrix.SE5b.20, Bp.Sp.5b,  
       SpBio20.matrix5b+SpBio.matrix.SE5b.20, Bp.Sp.5b,
       lwd=1.5, angle=90, code=3, length=0.03)

#Plotting the rank of the top 20 estimated in log-axis scale
LogBp.Sp.5b<-barplot(SpBio20.matrix5b, beside=T, horiz=T, log="x",
                     legend=rownames(SpBio20.matrix5b), xlim=c(0.1, 700),
                     args.legend=list(x="bottomright", cex=0.8, text.width=0.3, bty="n"),
                     col= brewer.pal(n=6, name="Set1"),font=3, xaxt="n",
                     xlab="Estimated detection per 100 km/year", ylab="")
title(ylab="Species", line=11, cex.lab=1)
marks<-c(0.1, 0.5, 1, 5, 10, 50, 100, 500)
axis(1,at=marks,labels=marks)

segments(SpBio20.matrix5b-SpBio.matrix.SE5b.20, LogBp.Sp.5b,
         SpBio20.matrix5b+SpBio.matrix.SE5b.20, LogBp.Sp.5b,
         lwd=1.5)

arrows(SpBio20.matrix5b-SpBio.matrix.SE5b.20, LogBp.Sp.5b,
       SpBio20.matrix5b+SpBio.matrix.SE5b.20, LogBp.Sp.5b,
       lwd=1.5, angle=90, code=3, length=0.03)


#Ranking the endangered species####
# Revaluing the name of species of conservation concern###
SpBio.df5b$Sp.CC<-revalue(SpBio.df5b$Sp,c("Tangara_peruviana" = "Tangara peruviana*",
                                          "Priodontes_maximus" = "Priodontes maximus*",
                                          "Leopardus_tigrinus" = "Leopardus tigrinus****",
                                          "Blastocerus_dichotomus" = "Blastocerus dichotomus*", 
                                          "Callithrix_aurita" = "Callithrix aurita****",
                                          "Tapirus_terrestris" = "Tapirus terrestris*",
                                          "Ramphastos_vitellinus" = "Ramphasto vitellinus***",
                                          "Tayassu_pecari" = "Tayassu pecari*",
                                          "Myrmecophaga_tridactyla" = "Myrmecophaga tridactyla*",
                                          "Puma_yagouaroundi" =      "Puma yagouaroundi**",
                                          "Chrysocyon_brachyurus" = "Chrysocyon brachyurus**",
                                          "Leopardus_wiedii" = "Leopardus wiedii**",
                                          "Alouatta_guariba_clamitans" = "Alouatta guariba clamitans**",
                                          "Leopardus_colocolo" = "Leopardus colocolo**", 
                                          "Leopardus_geoffroyi" = "Leopardus geoffroyi**", 
                                          "Lycalopex_vetulus" = "Lycalopex vetulus**",
                                          "Ozotoceros_bezoarticus" = "Ozotoceros bezoarticus**",
                                          "Panthera_onca" = "Panthera onca**",
                                          "Sapajus_cay" = "Sapajus cay**",
                                          "Puma_concolor" = "Puma concolor**"))

# Creating a matrix for the top 20 estimated species of conservation concern ####
CC.matrix5b<- tapply(SpBio.df5b$Estimate40, list(SpBio.df5b$Bio, SpBio.df5b$Sp.CC), mean)
CC.matrix5b<-CC.matrix5b[, c("Tangara peruviana*","Priodontes maximus*",
                             "Leopardus tigrinus****","Blastocerus dichotomus*", 
                             "Callithrix aurita****","Tapirus terrestris*",
                             "Ramphasto vitellinus***","Tayassu pecari*",
                             "Myrmecophaga tridactyla*", "Puma yagouaroundi**",
                             "Chrysocyon brachyurus**", "Leopardus wiedii**",
                             "Alouatta guariba clamitans**","Leopardus colocolo**", 
                             "Leopardus geoffroyi**",  "Lycalopex vetulus**",
                             "Ozotoceros bezoarticus**","Panthera onca**",
                             "Sapajus cay**", "Puma concolor**")]
CC.matrix.colmax5b<- apply(CC.matrix5b, MARGIN=2, FUN=max, na.rm= T)
CC.matrix5b<- CC.matrix5b[, order(CC.matrix.colmax5b)]

# Creating a matrix for the standard error of the top 20 estimated species of conservation concern####
CC.matrix.SE5b<- tapply(SpBio.df5b$SE40, list(SpBio.df5b$Bio, SpBio.df5b$Sp.CC), mean)
CC.matrix.SE5b<-CC.matrix.SE5b[, c("Tangara peruviana*","Priodontes maximus*",
                                   "Leopardus tigrinus****","Blastocerus dichotomus*", 
                                   "Callithrix aurita****","Tapirus terrestris*",
                                   "Ramphasto vitellinus***","Tayassu pecari*",
                                   "Myrmecophaga tridactyla*", "Puma yagouaroundi**",
                                   "Chrysocyon brachyurus**", "Leopardus wiedii**",
                                   "Alouatta guariba clamitans**","Leopardus colocolo**", 
                                   "Leopardus geoffroyi**",  "Lycalopex vetulus**",
                                   "Ozotoceros bezoarticus**","Panthera onca**",
                                   "Sapajus cay**", "Puma concolor**")]
CC.matrix.SE5b<-CC.matrix.SE5b[, order(CC.matrix.colmax5b)]

# Plotting the rank of species of conservation concern
par(mfrow=c(1,1), mar= c(4.1, 12, 4.1, 4.1), las=1)
Bp.Sp.CC.5b<-barplot(CC.matrix5b, beside=T, horiz=T, font=3, xlim=c(0,25),
                     legend=rownames(CC.matrix5b), 
                     args.legend=list(x="bottomright", cex=0.9, text.width=3, bty="n"), 
                     col= brewer.pal(n=6, name="Set1"),
                     xlab="Total of detected carcasses",
                     ylab="")
title(ylab="Species", line=11, cex.lab=1)

segments(CC.matrix5b - CC.matrix.SE5b, Bp.Sp.CC.5b,
         CC.matrix5b + CC.matrix.SE5b, Bp.Sp.CC.5b,
         lwd=1.5)

arrows(CC.matrix5b - CC.matrix.SE5b, Bp.Sp.CC.5b,
       CC.matrix5b + CC.matrix.SE5b, Bp.Sp.CC.5b,
       lwd=1.5, angle=90, code=3, length=0.03)

#Ranking the species of conservation concern per biome####
par(mfrow=c(2,2), mar=c(4.1, 12.1, 1.1, 3.1), las=1)

#Species of conservation concern in the Atlantic forest####
AFdata.CC5b<-SpBio.df5b[SpBio.df5b$Bio=="Atlantic_forest",]

# Creating a matrix for the top 20 estimated species of conservation concern in the Atlantic forest####
AFCC.matrix5b<- tapply(AFdata.CC5b$Estimate40, list(AFdata.CC5b$Bio, AFdata.CC5b$Sp.20), mean)
AFCC.matrix5b<- AFCC.matrix5b[, c("Leopardus colocolo","Tangara peruviana",
                                  "Leopardus wiedii", "Myrmecophaga tridactyla",
                                  "Leopardus tigrinus","Callithrix aurita",
                                  "Ramphastos vitellinus","Puma yagouaroundi",
                                  "Chrysocyon brachyurus")]
AFCC.matrix5b.col<-apply(AFCC.matrix5b, MARGIN=2, FUN=max, na.rm= T)
AFCC.matrix5b<- AFCC.matrix5b[, order(AFCC.matrix5b.col)]

# Creating a matrix for the standard error of the top 20 estimated species of conservation concern in the Atlantic forest####
AFCC.matrix5b.SE<- tapply(AFdata.CC5b$SE40, list(AFdata.CC5b$Bio, AFdata.CC5b$Sp.20), mean)
AFCC.matrix5b.SE<-AFCC.matrix5b.SE[, c("Leopardus colocolo","Tangara peruviana",
                                       "Leopardus wiedii", "Myrmecophaga tridactyla",
                                       "Leopardus tigrinus","Callithrix aurita",
                                       "Ramphastos vitellinus","Puma yagouaroundi",
                                       "Chrysocyon brachyurus")]
AFCC.matrix5b.SE<-AFCC.matrix5b.SE[, order(AFCC.matrix5b.col)]

# Plotting the rank of species of conservation concern in the Atlantic forest
Bp.AFCC.5b<-barplot(AFCC.matrix5b, col="#0066CC",  
                    xlim=c(0, 3), beside=T, horiz=T,
                    xlab="Estimated detection per 100 km/year", xaxt="n", font=3,
                    ylab="")
title(ylab="Species", line=11, cex.lab=1)
marks<-c(0, 0.5, 1, 1.5, 2, 2.5, 3)
axis(1,at=marks,labels=marks)

segments(AFCC.matrix5b-AFCC.matrix5b.SE, Bp.AFCC.5b,  
         AFCC.matrix5b+AFCC.matrix5b.SE, Bp.AFCC.5b,
         lwd=1.5)

arrows(AFCC.matrix5b-AFCC.matrix5b.SE, Bp.AFCC.5b,  
       AFCC.matrix5b+AFCC.matrix5b.SE, Bp.AFCC.5b,
       lwd=1.5, angle=90, code=3, length=0.03)

#Species of conservation concern in Cerrado####
CEdata.CC5b<-SpBio.df5b[SpBio.df5b$Bio=="Cerrado",]

# Creating a matrix for the top 20 estimated species of conservation concern in Cerrado####
CECC.matrix5b<- tapply(CEdata.CC5b$Estimate40, list(CEdata.CC5b$Bio, CEdata.CC5b$Sp.20), mean)
CECC.matrix5b<- CECC.matrix5b[, c("Priodontes maximus","Leopardus tigrinus",
                                  "Blastocerus dichotomus","Sapajus cay",
                                  "Tapirus terrestris","Puma yagouaroundi",
                                  "Tayassu pecari","Puma concolor", 
                                  "Ozotoceros bezoarticus","Chrysocyon brachyurus",
                                  "Lycalopex vetulus","Myrmecophaga tridactyla")]
CECC.matrix5b.col<-apply(CECC.matrix5b, MARGIN=2, FUN=max, na.rm= T)
CECC.matrix5b<- CECC.matrix5b[, order(CECC.matrix5b.col)]

# Creating a matrix for the standard error of the top 20 estimated species of conservation concern in Cerrado####
CECC.matrix5b.SE<- tapply(CEdata.CC5b$SE40, list(CEdata.CC5b$Bio, CEdata.CC5b$Sp.20), mean)
CECC.matrix5b.SE<-CECC.matrix5b.SE[, c("Priodontes maximus","Leopardus tigrinus",
                                       "Blastocerus dichotomus","Sapajus cay",
                                       "Tapirus terrestris","Puma yagouaroundi",
                                       "Tayassu pecari","Puma concolor", 
                                       "Ozotoceros bezoarticus","Chrysocyon brachyurus",
                                       "Lycalopex vetulus","Myrmecophaga tridactyla")]
CECC.matrix5b.SE<-CECC.matrix5b.SE[, order(CECC.matrix5b.col)]

# Plotting the rank of species of conservation concern in Cerrado
Bp.CECC.5b<-barplot(CECC.matrix5b,col="#9900CC",  
                    xlim=c(0, 25), beside=T, horiz=T,
                    xlab="Estimated detection per 100 km/year", font=3,
                    ylab="")
title(ylab="Species", line=11, cex.lab=1)

segments(CECC.matrix5b-CECC.matrix5b.SE, Bp.CECC.5b,  
         CECC.matrix5b+CECC.matrix5b.SE, Bp.CECC.5b,
         lwd=1.5)

arrows(CECC.matrix5b-CECC.matrix5b.SE, Bp.CECC.5b,  
       CECC.matrix5b+CECC.matrix5b.SE, Bp.CECC.5b,
       lwd=1.5, angle=90, code=3, length=0.03) 

#Species of conservation concern in Pampa####
PPdata.CC5b<-SpBio.df5b[SpBio.df5b$Bio=="Pampa",]

# Creating a matrix for the top 20 estimated species of conservation concern in Pampa####
PPCC.matrix5b<- tapply(PPdata.CC5b$Estimate40, list(PPdata.CC5b$Bio, PPdata.CC5b$Sp.20), mean)
PPCC.matrix5b<- PPCC.matrix5b[, c("Alouatta guariba_clamitans","Leopardus wiedii",
                                  "Leopardus geoffroyi")]
PPCC.matrix5b.col<-apply(PPCC.matrix5b, MARGIN=2, FUN=max, na.rm= T)
PPCC.matrix5b<- PPCC.matrix5b[, order(PPCC.matrix5b.col)]

# Creating a matrix for the standard error of the top 20 estimated species of conservation concern in Pampa####
PPCC.matrix5b.SE<- tapply(PPdata.CC5b$SE40, list(PPdata.CC5b$Bio, PPdata.CC5b$Sp.20), mean)
PPCC.matrix5b.SE<-PPCC.matrix5b.SE[,c("Alouatta guariba_clamitans","Leopardus wiedii",
                                      "Leopardus geoffroyi")]
PPCC.matrix5b.SE<-PPCC.matrix5b.SE[, order(PPCC.matrix5b.col)]

# Plotting the rank of species of conservation concern in Pampa
Bp.PPCC.5b<-barplot(PPCC.matrix5b, col="#FF9933", #brewer.pal(n=6, name="Set1"), 
                    xlim=c(0, 20), beside=T, horiz=T,
                    xlab="Estimated detection per 100 km/year", font=3,
                    ylab="")
title(ylab="Species", line=11, cex.lab=1)

segments(PPCC.matrix5b-PPCC.matrix5b.SE, Bp.PPCC.5b,  
         PPCC.matrix5b+PPCC.matrix5b.SE, Bp.PPCC.5b,
         lwd=1.5)

arrows(PPCC.matrix5b-PPCC.matrix5b.SE, Bp.PPCC.5b,  
       PPCC.matrix5b+PPCC.matrix5b.SE, Bp.PPCC.5b,
       lwd=1.5, angle=90, code=3, length=0.03)

#Species of conservation concern in Pantanal####
PTdata.CC5b<-SpBio.df5b[SpBio.df5b$Bio=="Pantanal",]

# Creating a matrix for the top 20 estimated species of conservation concern in Pantanal####
PTCC.matrix5b<- tapply(PTdata.CC5b$Estimate40, list(PTdata.CC5b$Bio, PTdata.CC5b$Sp.20), mean)
PTCC.matrix5b<- PTCC.matrix5b[, c("Blastocerus dichotomus","Tapirus terrestris",
                                  "Panthera onca","Tayassu pecari",
                                  "Puma yagouaroundi", "Myrmecophaga tridactyla")]
PTCC.matrix5b.col<-apply(PTCC.matrix5b, MARGIN=2, FUN=max, na.rm= T)
PTCC.matrix5b<- PTCC.matrix5b[, order(PTCC.matrix5b.col)]

# Creating a matrix for the standard error of the top 20 estimated species of conservation concern in Pantanal####
PTCC.matrix5b.SE<- tapply(PTdata.CC5b$SE40, list(PTdata.CC5b$Bio, PTdata.CC5b$Sp.20), mean)
PTCC.matrix5b.SE<-PTCC.matrix5b.SE[,c("Blastocerus dichotomus","Tapirus terrestris",
                                      "Panthera onca","Tayassu pecari",
                                      "Puma yagouaroundi", "Myrmecophaga tridactyla")]
PTCC.matrix5b.SE<-PTCC.matrix5b.SE[, order(PTCC.matrix5b.col)]

# Plotting the rank of species of conservation concern in Pantanal
Bp.PTCC.5b<-barplot(PTCC.matrix5b,col= "#FFFF33", #brewer.pal(n=1, name="Set1"), 
                    xlim=c(0, 10), beside=T, horiz=T,
                    xlab="Estimated detection per 100 km/year", font=3,
                    ylab="", 
                    legend= c("Pantanal", "Pampa", "Cerrado", "Atlantic forest"), #rownames(PTCC.matrix5b), 
                    args.legend=list(x="bottomright", cex=0.9, text.width=3, bty="n", fill=c("#0066CC", "#9900CC", "#FF9933", "#FFFF33")))
title(ylab="Species", line=11, cex.lab=1)

segments(PTCC.matrix5b-PTCC.matrix5b.SE, Bp.PTCC.5b,  
         PTCC.matrix5b+PTCC.matrix5b.SE, Bp.PTCC.5b,
         lwd=1.5)

arrows(PTCC.matrix5b-PTCC.matrix5b.SE, Bp.PTCC.5b,  
       PTCC.matrix5b+PTCC.matrix5b.SE, Bp.PTCC.5b,
       lwd=1.5, angle=90, code=3, length=0.03)


# Finding the number of recorded carcasses for a given species per biome ####
tapply(RK7Rp.final$Quantity, list(RK7Rp.final$Biome, RK7Rp.final$Species=="Didelphis_aurita"), sum)
tapply(RK7Rp.final$Quantity, list(RK7Rp.final$Biome, RK7Rp.final$Species=="Didelphis_albiventris"), sum)
tapply(RK7Rp.final$Quantity, list(RK7Rp.final$Biome, RK7Rp.final$Species=="Trachemys_dorbigni"), sum)

# Finding the number of estimated carcasses and the standard error for a given species per biome ####
tapply(SpBio.df5b$Estimate40, list(SpBio.df5b$Bio, SpBio.df5b$Sp=="Myrmecophaga_tridactyla"), sum)
tapply(SpBio.df5b$SE40, list(SpBio.df5b$Bio, SpBio.df5b$Sp=="Myrmecophaga_tridactyla"), sum)
tapply(SpBio.df5b$Estimate40, list(SpBio.df5b$Bio, SpBio.df5b$Sp=="Leopardus_geoffroyi"), sum)
tapply(SpBio.df5b$SE40, list(SpBio.df5b$Bio, SpBio.df5b$Sp=="Leopardus_geoffroyi"), sum)

#Calculating the estimates and SE per fauna group ####
E.group5b<-tapply(SpBio.df5b$Estimate40, SpBio.df5b$Fauna_group, sum, na.rm=T)
E.group5b
SE.group5b<-tapply(SpBio.df5b$SE40, SpBio.df5b$Fauna_group, sum, na.rm=T)
SE.group5b

#Ranking the top 15 species per Fauna group####
par(mfrow=c(1,3), mar= c(4.1, 12, 3.1, 4.1), las=1, cex=0.8)
#pdf("Rank_SpBio.pdf", paper="a4")
#Birds####
Bird.5b<-SpBio.df5b[SpBio.df5b$Fauna_group=="A",]

Bird.matrix5b<- tapply(Bird.5b$Estimate40, list(Bird.5b$Bio, Bird.5b$Sp.20), mean)
Bird.matrix5b.col<-apply(Bird.matrix5b, MARGIN=2, FUN=max, na.rm= T)
Bird.matrix5b<- Bird.matrix5b[, order(Bird.matrix5b.col)]
Bird.matrix5b<-Bird.matrix5b[,285:299]
Bird.matrix5b.SE<- tapply(Bird.5b$SE40, list(Bird.5b$Bio, Bird.5b$Sp.20), mean)
Bird.matrix5b.SE<-Bird.matrix5b.SE[, order(Bird.matrix5b.col)]
Bird.matrix5b.SE<-Bird.matrix5b.SE[, 285:299]

Bird.CA.5b<-barplot(Bird.matrix5b, beside=T, horiz=T, font=3, log="x",
                    col= brewer.pal(n=6, name="Set1"), xlim=c(0.1, 1000), 
                    xlab="", main= "a) Birds",
                    ylab="", xaxt="n")
title(ylab="Species", line=10, cex.lab=1.2)
marks<-c(0.1, 1, 10, 100, 1000)
axis(1,at=marks,labels=marks)

segments(Bird.matrix5b-Bird.matrix5b.SE, Bird.CA.5b,  
         Bird.matrix5b+Bird.matrix5b.SE, Bird.CA.5b,
         lwd=1.5)

arrows(Bird.matrix5b-Bird.matrix5b.SE, Bird.CA.5b,  
       Bird.matrix5b+Bird.matrix5b.SE, Bird.CA.5b,
       lwd=1.5, angle=90, code=3, length=0.03)

#Mammals#####
Mam.5b<-SpBio.df5b[SpBio.df5b$Fauna_group=="M",]

Mam.matrix5b<- tapply(Mam.5b$Estimate40, list(Mam.5b$Bio, Mam.5b$Sp.20), mean)
Mam.matrix5b.col<-apply(Mam.matrix5b, MARGIN=2, FUN=max, na.rm= T)
Mam.matrix5b<- Mam.matrix5b[, order(Mam.matrix5b.col)]
Mam.matrix5b<-Mam.matrix5b[,94:108]
Mam.matrix5b.SE<- tapply(Mam.5b$SE40, list(Mam.5b$Bio, Mam.5b$Sp.20), mean)
Mam.matrix5b.SE<-Mam.matrix5b.SE[, order(Mam.matrix5b.col)]
Mam.matrix5b.SE<-Mam.matrix5b.SE[, 94:108]

Mam.CA.5b<-barplot(Mam.matrix5b, beside=T, horiz=T, font=3, log="x",
                   col= brewer.pal(n=6, name="Set1"), xlim=c(0.1, 1000), 
                   xlab="Estimated detection per 100 km/year ", 
                   ylab="", xaxt="n", main= "b) Mammals")
#title(ylab="Species", line=11, cex.lab=1)
marks<-c(0.1, 1, 10, 100, 1000)
axis(1,at=marks,labels=marks)

segments(Mam.matrix5b-Mam.matrix5b.SE, Mam.CA.5b,  
         Mam.matrix5b+Mam.matrix5b.SE, Mam.CA.5b,
         lwd=1.5)

arrows(Mam.matrix5b-Mam.matrix5b.SE, Mam.CA.5b,  
       Mam.matrix5b+Mam.matrix5b.SE, Mam.CA.5b,
       lwd=1.5, angle=90, code=3, length=0.03)


#Reptiles#####
#par(oma=c(0,0,0,4)) # set margins par(oma = c(b, l, t, r))
Rep.5b<-SpBio.df5b[SpBio.df5b$Fauna_group=="H",]

Rep.matrix5b<- tapply(Rep.5b$Estimate40, list(Rep.5b$Bio, Rep.5b$Sp.20), mean)
Rep.matrix5b.col<-apply(Rep.matrix5b, MARGIN=2, FUN=max, na.rm= T)
Rep.matrix5b<- Rep.matrix5b[, order(Rep.matrix5b.col)]
Rep.matrix5b<-Rep.matrix5b[,37:51]
Rep.matrix5b.SE<- tapply(Rep.5b$SE40, list(Rep.5b$Bio, Rep.5b$Sp.20), mean)
Rep.matrix5b.SE<-Rep.matrix5b.SE[, order(Rep.matrix5b.col)]
Rep.matrix5b.SE<-Rep.matrix5b.SE[, 37:51]

#par(oma = c(4.1,12,3.1,10), new = TRUE)
Rep.CA.5b<-barplot(Rep.matrix5b, beside=T, horiz=T, font=3, log="x",
                   col= brewer.pal(n=6, name="Set1"), xlim=c(0.1, 1000), 
                   xlab="", main= "c) Reptiles",
                   ylab="", xaxt="n",
                   legend.text=rownames(Rep.matrix5b),
                   args.legend=list(x="bottomright", bty="n",
                                    inset=c(-0.50,0), xpd=T))
marks<-c(0.1, 1, 10, 100, 1000)
axis(1,at=marks,labels=marks)

segments(Rep.matrix5b-Rep.matrix5b.SE, Rep.CA.5b,  
         Rep.matrix5b+Rep.matrix5b.SE, Rep.CA.5b,
         lwd=1.5)

arrows(Rep.matrix5b-Rep.matrix5b.SE, Rep.CA.5b,  
       Rep.matrix5b+Rep.matrix5b.SE, Rep.CA.5b,
       lwd=1.5, angle=90, code=3, length=0.03)
#dev.off()

# Obtaining the fitted values on the original scale ####

# Plotting the relationship between the estimated detection rates and the LogBM ####
par(mfrow=c(1,1), mar= c(4.1, 4.1, 4.1, 4.1), las=1)
plot(SpBio.df5b$coef~SpBio.df5b$LogBM, xlab="Body mass (log)", ylab="Estimated detection",
     pch=16)

###############
# identifying typical small / medium / large species abundances ####
SpBio.df5b[identify(SpBio.df5b$coef~SpBio.df5b$LogBM, n=3),]

#Using 10 quantile, mean and 90 quantile to identify small, medium and large species####
summary(SpBio.df5b$adult_BM_g)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#3.40     26.36    120.85   3012.87   1033.74 225000.00

#Using the species in the 1st quartile as a small size species
#quantile(SpBio.df5b$adult_BM_g, 0.1) # 12.08  - no species at this size, but there is at 12
SpBio.df5b[SpBio.df5b$adult_BM_g==12 ,]
#                   coef         code                       Sp             Bio    LogBM IUCN BRL Fauna_group
# s(SpBio).279 -1.899147 s(SpBio).279 Hemithraupis_ruficapilla Atlantic_forest 2.484907   LC  LC           A
# s(SpBio).558  1.429273 s(SpBio).558     Sporophila_americana          Amazon 2.484907   LC  LC           A
#              TotalCarcasses Estimate40                                    SpBio      SE40    minSE40    maxSE40
# s(SpBio).279              1  0.2307836 Hemithraupis_ruficapilla_Atlantic_forest 0.1682189 0.06256466  0.3990025
# s(SpBio).558              3  6.4375315              Sporophila_americana_Amazon 5.6568680 0.78066354 12.0943995
#                    PSE Estimate50  Estimate60 adult_BM_g
# s(SpBio).279 72.89034 0.02086346 0.001886114         12
# s(SpBio).558 87.87325 0.58197031 0.052611695         12

RK.S5b<-data.frame(Car_velocity=seq(from=40, to=60, length=20), 
                   LogBM=2.484907, 
                   SpBio="Sporophila_americana_Amazon",
                   Gap_Size="<4",
                   LogSection=log(100*365),
                   S.Gs4plus="dummyLevel"
) 

#Using the species with the mean weight as a medium size species
summary(SpBio.df5b$adult_BM_g)[4] # 3012.871  - There is no species with that BM but there is at 3113
SpBio.df5b[SpBio.df5b$adult_BM_g==3113 ,]
#                  coef        code               Sp   Bio    LogBM IUCN BRL Fauna_group TotalCarcasses Estimate40
# s(SpBio).46 1.567064 s(SpBio).46 Boiruna_maculata Pampa 8.043342   NI  NI           H              7   12.45547
#                              SpBio     SE40  minSE40  maxSE40      PSE Estimate50 Estimate60 adult_BM_g
# s(SpBio).46 Boiruna_maculata_Pampa 7.089102 5.366369 19.54457 56.91557   9.973221   7.985658       3113

RK.M5b<-data.frame(Car_velocity=seq(from=40, to=60, length=20), 
                   LogBM=8.043342,
                   SpBio="Boiruna_maculata_Pampa",
                   Gap_Size="<4",
                   LogSection=log(100*365),
                   S.Gs4plus="dummyLevel"
)  

#Using the species in the 95% quantile as a large size species
SpBio.df5b[SpBio.df5b$adult_BM_g==quantile(SpBio.df5b$adult_BM_g, 0.95),] 
#                    coef         code                 Sp             Bio    LogBM IUCN BRL Fauna_group TotalCarcasses
# s(SpBio).302 -0.8842920 s(SpBio).302 Leopardus_pardalis Atlantic_forest 9.259131   LC  NI           M             11
# s(SpBio).303  0.3033384 s(SpBio).303 Leopardus_pardalis        Caatinga 9.259131   LC  NI           M              1
# s(SpBio).304 -0.4788869 s(SpBio).304 Leopardus_pardalis         Cerrado 9.259131   LC  NI           M              9
# s(SpBio).305  0.3748982 s(SpBio).305 Leopardus_pardalis        Pantanal 9.259131   LC  NI           M             11
#              Estimate40                              SpBio      SE40   minSE40  maxSE40      PSE Estimate50
# s(SpBio).302   1.203254 Leopardus_pardalis_Atlantic_forest 0.3556818 0.8475721 1.558936 29.56000   1.552504
# s(SpBio).303   3.945832        Leopardus_pardalis_Caatinga 3.6478034 0.2980288 7.593636 92.44700   5.091128
# s(SpBio).304   1.804773         Leopardus_pardalis_Cerrado 0.6969351 1.1078376 2.501708 38.61623   2.328616
# s(SpBio).305   4.238544        Leopardus_pardalis_Pantanal 2.5984047 1.6401389 6.836948 61.30419   5.468801
#               Estimate60 adult_BM_g
# s(SpBio).302   2.003126      10500
# s(SpBio).303   6.568852      10500
# s(SpBio).304   3.004508      10500
# s(SpBio).305   7.056145      10500

RK.L5b<-data.frame(Car_velocity=seq(from=40, to=60, length=20), 
                   LogBM=9.259131  ,
                   SpBio="Leopardus_pardalis_Pantanal",
                   Gap_Size="<4",
                   LogSection=log(100*365),
                   S.Gs4plus="dummyLevel"
)

pred.1.5b<-predict(gam.NB.random5b, newdata=RK.S5b, type="link")
pred.2.5b<-predict(gam.NB.random5b, newdata=RK.M5b, type="link")
pred.3.5b<-predict(gam.NB.random5b, newdata=RK.L5b, type="link")

pred.1R.5b<-predict(gam.NB.random5b, newdata=RK.S5b, type="response")
pred.2R.5b<-predict(gam.NB.random5b, newdata=RK.M5b, type="response")
pred.3R.5b<-predict(gam.NB.random5b, newdata=RK.L5b, type="response")

par(mfrow=c(1,1), mar= c(4.1, 8, 4.1, 4.1), cex=1.2, bty="n")

matplot(RK.S5b$Car_velocity, cbind(pred.1.5b, pred.2.5b, pred.3.5b), lty= 1, type= "l", lwd= 2,
        xlab="Velocity of survey (km/h)", ylab="Estimated detection per 100 km/year"
) # log scale

labs<-c("small", "medium", "large")
cols<-c("black", "red", "green")
legend("bottomleft", labs, col=cols, lty=1, cex=0.9, lwd=2, text.width=1) # log scale

matplot(RK.S5b$Car_velocity, cbind(pred.1R.5b, pred.2R.5b, pred.3R.5b), lty= 1, type= "l", lwd= 2,
        ylim=c(0,14), xlab="Survey speed (km/h)", ylab="Estimated detection per 100 km/year"
) # response scale

# Improving the same graph above
par(mfrow=c(1,1), mar= c(4.1, 7.1, 1.1, 4.1), cex=1.2, bty="n")

matplot(RK.S5b$Car_velocity, cbind(pred.1R.5b, pred.2R.5b, pred.3R.5b), lty= 1, type= "l", 
        lwd= 1.5, cex.lab=0.7, cex.axis=0.6,
        ylim=c(0,14), xlab="Survey speed (km/h)", ylab="Estimated detection per 100 km/year"
) # response scale

legend("topright", labs, col=cols, lty=1, lwd=1.5, cex=0.6, text.width=4, bty="n") # response scale

# Calculating the standard error for the estimates in the response scale ####
pred.1R.5b.se<-predict(gam.NB.random5b, newdata=RK.S5b, se.fit=T, type="response")
pred.1R.5b.se
pred.2R.5b.se<-predict(gam.NB.random5b, newdata=RK.M5b, se.fit=T, type="response")
pred.2R.5b.se
pred.3R.5b.se<-predict(gam.NB.random5b, newdata=RK.L5b, se.fit=T, type="response")
pred.3R.5b.se

#How many species with less than 0.3 carcasses estimated?
L3<-SpBio.df5b[SpBio.df5b$Estimate40<=0.3,]
str(L3)
summary(L3)
head(L3)
tail(L3)
#tapply(L3$Estimate40, L3$Fauna_group, mean, na.rm=T)

#Species Trachemys - carcass encounter per month per year
Trachemys<-RK7Rp.final[RK7Rp.final$Species=="Trachemys_dorbigni",]
Trachemys$NewDate<-as.Date(Trachemys$Date, "%d/%m/%Y")
Trachemys$Juliandate<-julian(as.Date(Trachemys$NewDate,"%d/%m/%Y"))
df<-data.frame(date = Trachemys$NewDate,
               year = as.numeric(format(Trachemys$NewDate, format = "%Y")),
               month = as.numeric(format(Trachemys$NewDate, format = "%m")),
               day = as.numeric(format(Trachemys$NewDate, format = "%d")))
Trachemys$Year<-df$year
Trachemys$Month<-df$month
Trachemys$Day<-df$day

mean.MY<- tapply(Trachemys$Quantity, list(Trachemys$Month, Trachemys$Year), mean, na.rm=T)
# matplot draws one line per column (year)
par(mfrow=c(1,1),mar= c(4.1, 4.1, 3.1, 3.1), las=1)
matplot(mean.MY, type= "l",bty= "n",
        ylim= c(0, 12), cex=1, xlab= "month", ylab= "Number of carcasses")

legend(x= "topleft", legend= colnames(mean.MY),
       bty= "n", # no bounding box for the legend
       col= 1:ncol(mean.MY),
       lty= 2:ncol(mean.MY),
       cex= 1,
       title= "Year")


