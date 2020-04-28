require(lme4)
require(lmerTest)
require(ggplot2)
require(doBy)
require(MuMIn)
require(broom)
require(data.table)
require(sjPlot)
require(e1071)
require(insight)
require(dplyr)
require(lattice)
require(Hmisc)
#=================================================================================================#
###READ IN DATA
grow1 <- read.csv(file="master6-15_grw_climNA.csv", sep=",",head=TRUE, na.string="na")
#OMIT NAs
grw <- na.omit(grow1)

#CREATE LOG-TRANSFORMED GROWTH DATA
skewness(grw$grwvol_m)
grw_log = log1p(grw$grwvol_m)
grw <- cbind(grw,grw_log)


#SUM ACROSS YEARS (2011 to 2012)
#grw1yr <- summaryBy(grwvol_m+grw_log ~ pop + ssp + ploidy + family + garden+MAT+MWMT+MCMT+
#                      TD+MAP+MSP+AHM+SHM+DD0+DD5+DDL18+DDG18+NFFD+bFFP+eFFP+FFP+PAS+EMT+EXT+
#                      MAR+Eref+CMD+RH+PPT_wt+PPT_sp+PPT_sm+PPT_at+DD0_wt+DD5_sm,data = grw,FUN = c(sum))
#setnames(grw1yr, old = c('grwvol_m.sum','grw_log.sum'),new = c('grwvol_m','grw_log'))
#SUMMARIZE 
exp_sum1 <- summaryBy(pop+grwvol_m~pop+ssp+ploidy+garden, data = grw,FUN=c(length,mean,sd))
exp_sum2 <- summaryBy(grwvol_m~ssp+garden, data = grw,FUN=c(mean,min,max,sd))
exp_sum3 <- summaryBy(pop+grwvol_m~pop+ssp+ploidy+garden+year, data = grw,FUN=c(length,mean,sd))

#APPENDIX S2 = exp_sum2
#write.csv(exp_sum1,file = "grw_exp_summary1.csv")
#write.csv(exp_sum2,file = "grw_exp_summary2.csv")
#=================================================================================================#

#Boxplot
rr <- ggplot(grw, aes(ssp,grwvol_m))
rr+geom_boxplot()+ facet_grid(.~garden)+ylim(-0.5,10)+xlab("subspecies")+ylab("growth")+
  theme_linedraw()

#Barplot (Fig.1)
ss <- ggplot(exp_sum3, aes(x=ssp,y=grwvol_m.mean,fill=factor(year, levels=c("y12","y11" ))))
ss + geom_bar(stat = "summary", fun.y = "mean")+facet_grid(.~garden)+ xlab("subspecies")+
  ylab(bquote('growth'~(m^3)))+theme_bw(base_size = 15)+guides(fill=guide_legend(title=NULL))


#=================================================================================================#
##SUMMARY BY POPULATION FOR ASSESSING CORRELATION WITH CLIMATE VARS
popgrw <- summaryBy(grw_log+MAT+MWMT+MCMT+TD+MAP+MSP+AHM+SHM+DD0+DD5+DDL18+DDG18+NFFD+bFFP+eFFP+FFP+
                      PAS+EMT+EXT+MAR+Eref+CMD+RH+PPT_wt+PPT_sp+PPT_sm+PPT_at+DD0_wt+DD5_sm ~ pop+garden+ssp, FUN = c(mean), data=grw)
setnames(popgrw, old = c('grw_log.mean','MAT.mean','MWMT.mean','MCMT.mean','TD.mean','MAP.mean','MSP.mean','AHM.mean','SHM.mean','DD0.mean','DD5.mean','DDL18.mean','DDG18.mean','NFFD.mean','bFFP.mean','eFFP.mean','FFP.mean','PAS.mean','EMT.mean','EXT.mean','MAR.mean','Eref.mean','CMD.mean','RH.mean','PPT_wt.mean','PPT_sp.mean','PPT_sm.mean','PPT_at.mean','DD0_wt.mean','DD5_sm.mean'),
         new = c('grw_log','MAT','MWMT','MCMT','TD','MAP','MSP','AHM','SHM','DD0','DD5','DDL18','DDG18','NFFD','bFFP','eFFP','FFP','PAS','EMT','EXT','MAR','Eref','CMD','RH','PPT_wt','PPT_sp','PPT_sm','PPT_at','DD0_wt','DD5_sm'))

popgrw_sum <- summaryBy(grw_log+MAT+MWMT+MCMT+TD+MAP+MSP+AHM+SHM+DD0+DD5+DDL18+DDG18+NFFD+bFFP+eFFP+FFP+
                          PAS+EMT+EXT+MAR+Eref+CMD+RH+PPT_wt+PPT_sp+PPT_sm+PPT_at+DD0_wt+DD5_sm ~ pop+garden+ssp, FUN = c(mean), data=grw)

#write.csv(x = popgrw_sum,file = "data_pop.csv")

data_pop=with(popgrw, data.frame(MAT,MWMT,MCMT,TD,MAP,MSP,AHM,SHM,DD0,DD5,DDL18,DDG18,NFFD,bFFP,eFFP,FFP,PAS,EMT,EXT,MAR,Eref,CMD,RH,PPT_wt,PPT_sp,PPT_sm,PPT_at,DD0_wt,DD5_sm,grw_log)) 
cor <- cor(data_pop)

####VARS CHOSEN r > 0.3 = MCMT, DD0, DDL18, PAS, MAT
##MCMT AND PAS HIGHEST NON COLINEAR VARS

#=================================================================================================#
###LMER
grw1 <- lmer(grw_log ~ ssp + PAS + (1|garden:year) + (1|pop:(garden:year)),REML=TRUE, data=grw)
grw2 <- lmer(grw_log ~ ssp + MCMT + (1|garden:year) + (1|pop:(garden:year)),REML=TRUE , data=grw)
grw3 <- lmer(grw_log ~ ssp + MCMT + PPT_sm + (1|garden:year) + (1|pop:(garden:year)),REML=TRUE , data=grw)
grw_wo_clim <- lmer(grw_log ~ ssp + (1|garden:year) + (1|pop:(garden:year)),REML=TRUE, data=grw)

##MODEL SELECTION WITH LRT
anova(grw1,grw2,grw3,grw_wo_clim)

### RESULTS
summary(grw1)
ranova(grw1)
r.squaredGLMM(grw1)
r.squaredGLMM(grw_wo_clim)
grw1
#Change to lmerMod for tables (Table 2)
class(grw1) <- "lmerMod"
tab_model(grw1,p.val = "kr", show.icc = FALSE)

##LMER OUTPUT and MAKE DATAFRAME
grw2_output <- augment(grw1)
grw2_output <- merge(grw2_output,grw)
grw2_summary <- summaryBy(grw_log+grwvol_m+.fitted+.fixed+PAS~
                            pop+ssp+garden+year,data=grw2_output,FUN = c(mean))
setnames(grw2_summary, old = c('grw_log.mean','grwvol_m.mean','.fitted.mean','.fixed.mean','PAS.mean'),
         new = c('grw_log','grwvol_m','.fitted','.fixed','PAS'))
#write.csv(grw2_summary,file = "grw_lmer_output.csv")

##MODEL FIT
#nn <- ggplot(grw2_summary, aes(grw_log,.fixed))
#nn + theme_bw()+stat_smooth(method=lm,se=FALSE,linetype=4, color="gray")+
#  geom_point(aes(shape=ssp))+labs(x="Observed", y="Predicted",check_overlap=TRUE)+
#  scale_shape_manual(values=c(20,23))
#cor.test(formula=~grw_log+.fixed,data = grw2_summary)

###EXAMINE RESIDUALS
residuals <- resid(grw1) 
summary(residuals)
hist(residuals)
plot(residuals)


###SEED-GROWTH COMPARISON
##Fig.3 (observed grw vs observed seed yield)
comp <- read.csv(file="comp.csv", sep=",", head=TRUE)
ll<-ggplot(comp,aes(x=grwvol_m,y=weight))+geom_point(aes(shape = ssp),size = 2)+
  stat_smooth(method=lm,se=FALSE,linetype=4,color="gray")+theme_bw()+
  xlab(bquote('growth'~(m^3)))+ylab("seed yield (g)")+scale_shape(solid = FALSE)+
  scale_x_continuous(breaks = c(0,0.5,1,1.5,2,2.5,3,3.5)) + 
  scale_y_continuous(breaks = c(0,50,100,150,200,250,300,350,400)) + theme_bw(base_size = 15)
ll
cor.test(formula = ~ grwvol_m+weight, data = comp)

###GROWTH-SEED LINEAR MODEL
lmALL <- lm(weight~grwvol_m,data=comp)
summary(lmALL)
tab_model(lmALL)


#Relationship by garden
#mm<-ggplot(comp,aes(x=grwvol_m,y=weight))+geom_point(aes(shape = ssp))+
#  facet_grid(garden~.)+stat_smooth(method=lm,se=FALSE,linetype=4,color="gray",size=1)+theme_bw()+
#  xlab("growth")+ylab("seed yield")+scale_shape(solid = FALSE)+geom_text(aes(label=pop))
#mm

#GROWTH-SEED LINEAR MODEL
lmObserved <- lm(weight~grwvol_m + ssp + garden,data = comp)
lmObs_resid <- resid(lmObserved)

summary(lmObserved)
tab_model(lmObserved)

###EXTRACT RANDOM EFFECTS
re_pop <- ranef(grw1,condVar=TRUE, whichel = "pop:(garden:year)")
re_env <- ranef(grw1,condVar=TRUE, whichel = "garden:year")
#Fig.S1
dotplot(re_env)
#DELIST to EXAMINE POP INTERACTION
re_pop1 <- unlist(re_pop)
re_pop2 <- as.vector(re_pop1)

###EXAMINE GXE
#GROWTH (Fig.2a)
mm <- ggplot(exp_sum1, aes(y=grwvol_m.mean,x=garden,group=pop))
mm + geom_point()+ geom_line(alpha=0.4)+facet_grid(.~ssp) + 
  ylab(bquote('growth'~(m^3))) + theme_bw(base_size = 15)#+geom_text(aes(label=pop))

###SUMMARIZE FOR GxE
#BIND GxE RANDOM EFFECTS and SURVIVAL
grw2_summary_re <- cbind(grw2_summary,re_pop2)

#AVG YEAR
grw2_summary_re <- summaryBy(grw_log+grwvol_m+.fitted+.fixed+PAS
                          +re_pop2~pop+ssp+garden,data=grw2_summary_re,FUN = c(mean))

#Fig 4 Relationship between PAS and growth at each garden
ww <- ggplot(grw2_summary_re, aes(PAS.mean,grw_log.mean,color=ssp,shape=ssp))
ww+ geom_point()+stat_smooth(method=lm,se=FALSE,linetype=4,fill=NA)+theme_bw(base_size = 15)+
  xlab("Precipitation as Snow")+ylab("log(Growth)")+facet_grid(.~garden)

#SUMMARIZE TO GARDEN EFFECTS 
Ogrw_GxE <- dplyr::filter(grw2_summary_re, garden!="Ephraim",garden!="Majors")
Mgrw_GxE <- dplyr::filter(grw2_summary_re, garden!="Ephraim",garden!="Orchard")
Egrw_GxE <- dplyr::filter(grw2_summary_re, garden!="Majors",garden!="Orchard")

#GARDEN-SSP CORRELATION
with(subset(Ogrw_GxE, ssp==c("tridentata")),cor.test(PAS.mean, grw_log.mean))
with(subset(Ogrw_GxE, ssp==c("wyomingensis")),cor.test(PAS.mean, grw_log.mean))
with(subset(Mgrw_GxE, ssp==c("tridentata")),cor.test(PAS.mean, grw_log.mean))
with(subset(Mgrw_GxE, ssp==c("wyomingensis")),cor.test(PAS.mean, grw_log.mean))
with(subset(Egrw_GxE, ssp==c("tridentata")),cor.test(PAS.mean, grw_log.mean))
with(subset(Egrw_GxE, ssp==c("wyomingensis")),cor.test(PAS.mean, grw_log.mean))

#SEED YIELD POP X GARDEN INTERACTIONS (Fig.2b)
nn <- ggplot(comp, aes(y=weight,x=garden,group=pop))
nn + geom_point()+ geom_line(alpha=0.4)+facet_grid(.~ssp)+ 
  ylab("seed yield (g)") +theme_bw(base_size = 15) #+geom_text(aes(label=pop))

#GROWTH-MORTALITY 
#read in survival data, subset Ephraim
#surv <- read.csv(file = "surv_pop.csv",sep=",", head=TRUE)
#comp_surv <- read.csv(file = "comp_surv.csv",sep=",", head=TRUE)
#comp_surv <- cbind(comp_surv, surv$fitted, surv$a.hat,surv$observed)
#comp_cor <- comp_surv[c(5:20)]
#cor_traits <- cor(comp_cor)
