library(riskRegression)
library(ggfortify)
library(survminer)
library(plotly)
library(gridExtra)

data(Melanoma)
set.seed(234)
Melanoma$id=sample(1:205)
str(Melanoma)
table(Melanoma$event)

Melanoma$event=factor(Melanoma$event, levels = c("censored"
                                                 ,"death.malignant.melanoma",
                                                 "death.other.causes"),
                      labels = c("censurado","murio",
                                 "otra causa"))

Melanoma$sex=factor(Melanoma$sex, levels = c("Female","Male"),
                    labels = c("Mujer","Hombre"))



ggplotly(
  Melanoma %>%
    mutate(
      text = paste("id = ", id,
                   "<br>", "time= ", time,
                   "<br>", "status = ", status,
                   "<br>", " age= ", round(age, 2))
    ) %>%
    ggplot(aes(x = id, y = time, text = text)) +
    geom_linerange(aes(ymin = 0, ymax = time)) +
    geom_point(aes(shape = status, color = status), stroke = 1, cex = 2) +
    scale_shape_manual(values = c(1, 3, 4)) +
    labs(y = "Tiempo (días)", x = "id") + coord_flip() + theme_classic(),
  tooltip = "text"
)







###################################################
### 1.1
### Comparación no paramétrica de CIF


library(cmprsk)
cif<-cuminc(ftime = Melanoma$time, fstatus = Melanoma$status, group=Melanoma$sex)
plot(cif,col=1:4,xlab="Days")
cif






#########################################
#### 1.2 Cause-specific hazard regression
#### Regresión de riesgo de causa específica




### Regresión de Cox 

csh<-coxph(Surv(time,status==1)~sex+invasion,data=Melanoma)
summary(csh)





library("riskRegression")
library("prodlim")


CSH<-CSC(Hist(time,status)~sex+age+invasion,data=Melanoma)
CSH


library("pec")

pec::predictEventProb(CSH,cause=1,newdata=data.frame(age=40,invasion=factor("level.2",levels=levels(Melanoma$invasion)),
                                                     sex=factor("Hombre",levels=levels(Melanoma$sex))),
                      time=c(1000,2000,3000))




SH <- FGR(Hist(time,status)~sex+age+invasion,data=Melanoma)
SH





cov<-model.matrix(~sex+age+invasion,data=Melanoma)[,-1]
cov
crr.model<-crr(Melanoma$time,Melanoma$status,cov1=cov)
crr.model

### 1.4 Model prediction
### Predicción del modelo


newdata<-data.frame(sex=factor(c("Hombre","Hombre","Mujer"),
                               levels=levels(Melanoma$sex)),age=c(52,32,59),
                    invasion=factor(c("level.2","level.1","level.2"),
                                    levels=levels(Melanoma$invasion)))
newdata


dummy.new<-model.matrix(~sex+age+invasion,data=newdata)[,-1]
dummy.new


pred<-predict(crr.model,dummy.new)
plot(pred,lty=1:3,col=1:3,xlab="Días",ylab="Cumulative incidence function")
legend("topleft",c("Hombre,age=52,invasion2","Hombre,age=32,
                   invasion1","Female,age=59,invasion2"),lty=1:3,col=1:3)


#### Análisis de supervivencia en presencia de riesgos competitivos

reg<-riskRegression(Hist(time, status) ~ sex + age +invasion, data = Melanoma, cause = 1,link="prop")
reg
plot(reg,newdata=newdata)


### 1.5 Model diagnostic
### Modelo de diagnóstico


checkdata<-data.frame(sex=factor(c("Hombre","Hombre","Hombre"),
                                 levels=levels(Melanoma$sex)),
                      age=c(46,46,46),invasion=factor(c("level.0","level.1","level.2"),
                                                      levels=levels(Melanoma$invasion)))
checkdata
plot(reg,newdata=checkdata,lty=1:3,col=1:3)
text(2000,1,"Covariates sex='Hombre'; age=52")
legend("topleft",c("invasion.level0","invasion.level1","invasion.level2"),lty=1:3,col=1:3)




crr.time<-crr(Melanoma$time,Melanoma$status,cov1=cov,
              cov2=cov[,1],tf=function(t) t)
summary(crr.time)



# Regresión de riesgos competitivo


reg.time<-riskRegression(Hist(time, status) ~ sex + age +
                           strata(invasion), data = Melanoma, cause = 1,link="prop")

plotEffects(reg.time,formula=~invasion)


par(mfrow=c(2,2))
for(j in 1:ncol(crr.model$res)) {
  scatter.smooth(crr.model$uft, crr.model$res[,j],
                 main =names(crr.model$coef)[j],
                 xlab = "Failure time",
                 ylab ="Schoenfeld residuals")
}





crr.time2<-crr(Melanoma$time,Melanoma$status,cov1=cov,cov2=cov[,3],tf=function(t) t)
crr.time2

crr.time3<-crr(Melanoma$time,Melanoma$status,cov1=cov,cov2=cbind(cov[,3], cov[,3]),tf=function(t) cbind(t,t^2),)
crr.time3




















































