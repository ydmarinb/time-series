library(forecast)
library(TSA)
library(strucchange)
library(fANCOVA)
library(FitAR)
library(lmtest)
library(car)
Datos <- read.table(file.choose(),header=T,sep=";",skip=11,dec=",",colClasses=c(rep("NULL",11),"numeric",rep("NULL",6)))
Datos <- ts(Datos,freq=4,start=c(2000,1))

#--------------------------------Funciones---------------------------------

#Creando funci?on para extraer correctamente los valores deltai
factoresdeltai=function(descom,s,estacionini){ if(estacionini==1){ deltasi=descom$figure }
  if(estacionini!=1){ j=estacionini;deltasi=c(descom$figure[(s-j+2):s],descom$figure[1:(s-j+1)]) } 
  deltasi
}

#Creando funcion usuario crit.inf.resid() para calcular C_n?*(p)
crit.inf.resid=function(residuales,n.par,AIC="TRUE")
{ if(AIC=="TRUE"){
  #Calcula AIC
  CI=log(mean(residuales^2))+2*n.par/length(residuales)
} 
  if(AIC=="FALSE"){ 
    #Calcula BIC 
    CI=log(mean(residuales^2))+n.par*log(length(residuales))/length(residuales) 
  }
  CI
}


#Funci?on para calcular la amplitud de los I.P 
amplitud=function(LIP,LSP){ 
  a=LSP-LIP
  am=mean(a)
  am 
}

#Funci?on para calcular la cobertura de los I.P 
cobertura=function(real,LIP,LSP){
  I=ifelse(real>=LIP & real<=LSP,1,0) 
  p=mean(I)
  p 
}


BP.LB.test=function(serie,maxlag,type="Ljung"){ 
aux=floor(maxlag/6);
X.squared=c(rep(NA,aux))
 df=c(rep(NA,aux)) 
p.value=c(rep(NA,aux)) 
for(i in 1:aux){
test=Box.test(serie,lag=(6*i),type=type) 
X.squared[i]=test[[1]]
df[i]=test[[2]]
p.value[i]=test[[3]] 
}
 lag=6*c(1:aux)  
teste=as.data.frame(cbind(X.squared,df,p.value)) 
rownames(teste)=lag 
teste } 


pruebaDW1=function(modelo){
dwneg=durbinWatsonTest(modelo,max.lag=1,method="normal",alternative="negative")
dwpos=durbinWatsonTest(modelo,max.lag=1,method="normal",alternative="positive")
res=data.frame(1,dwneg$r,dwneg$dw,dwpos$p,dwneg$p)
names(res)=c("lag","rho estimado","Estadístico D-W","VP rho>0","VP rho<0")
res}



 
#------------------------------Grafica serie-------------------------------------------------------------------------------

win.graph()
g1 <- plot(Datos,main="Serie original",ylim=c(min(Datos),max(Datos)))
trans <- log(Datos)
win.graph()
g2 <- plot(trans,main="Logar?tmo  serie",ylim=c(min(trans),max(trans)))
acf(Datos,lag.max=m,ci.type="ma",ci.col=2,col=4)

#------------------Construcción de modelo------------------------------------------------
n <- length(Datos)-4
m <- round(n/4,0)
t <- 1:n
yt <- ts(Datos[t],freq=4,start=c(2000,1)) #valores de la serie en muestra de ajuste
lnyt <- log(yt)
trimestre <- seasonaldummy(yt) #Matriz con las 3 primeras variables Indicadoras 
#Separando una a una las 3 variables indicadoras
I1 <- trimestre[,1]
I2 <- trimestre[,2]
I3 <- trimestre[,3]
#Valores ajustados



#Modelo log-polinomial estacional
mod1 <- lm(lnyt~t+I(t^2)+I(t^3)+I(t^4)+I1+I2+I3)

yhatmod1=ts(exp(fitted(mod1))*exp(summary(mod1)$sigma^2/2),freq=4,start=c(2000,1))

summary(mod1)


#-----------------------------------Graficos modelo--------------------------------

#Grafica modelo 1

plot(Datos,main="Serie real y ajustada modelo1")
lines(yhatmod1,col=2,lwd=2)
legend("topleft",legend=c("Original","Ajustado modelo 1"),lty=1,col=c(1,2))

#------------------------------AIC y BIC---------------------------------------
et1 <- yt-yhatmod1
AIC=exp(crit.inf.resid(residuales=et1,n.par=8))
BIC=exp(crit.inf.resid(residuales=et1,n.par=8,AIC="FALSE"))
tablacriter=cbind(Log_polinomial=c(AIC,BIC))
rownames(tablacriter)=c("AIC","BIC")
tablacriter 
#-------------------------------------validacion cruzada----------------------

tnuevo <- (n+1):length(Datos) #?ndice de tiempo en los pron?sticos
ytnuevo=ts(Datos[tnuevo],freq=4,start=c(2016,4)) #valores reales de la serie en #muestra de validaci?on cruzada
trimestre.nuevo=seasonaldummy(yt,h=4) #matriz con las 11 
I1n <- trimestre.nuevo[,1]
I2n <- trimestre.nuevo[,2]
I3n <- trimestre.nuevo[,3]

#--------------------------------------------Pronosticos------------------------

pronmod1=exp(predict(mod1,newdata=data.frame(t=tnuevo,I1=I1n,I2=I2n,I3=I3n),interval="prediction",level=0.95))*exp(summary(mod1)$sigma^2/2)
pronmod1=ts(pronmod1,freq=4,start=c(2016,4))
pronmod1
accuracy(pronmod1[,1],ytnuevo)


plot(Datos)
lines(yhatmod1,col=2,lwd=2)
lines(pronmod,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 1","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))



#--------------------------------Graficos Residuales----------------------------------------

#residuales modelo 1
win.graph()
plot.ts(residuals(mod1),ylim=c(min(residuals(mod1),-2*summary(mod1)$sigma),
max(residuals(mod1),2*summary(mod1)$sigma)),
main="Residuos vs. t\nModelo Log polinomial grado cuatro con indicadoras")
abline(h=0,col=2)




plot(fitted(mod1),residuals(mod1),ylim=c(min(residuals(mod1),-2*summary(mod1)$sigma),
max(residuals(mod1),2*summary(mod1)$sigma)),
main="Residuos vs. ajustados\nModelo exponencial cuadr?atico estacional")
abline(h=0,col=2)
abline(h=c(-2*summary(mod1)$sigma,2*summary(mod1)$sigma),col=2)


#---------------------Postulacion de modelos con errores ARMA(p,q)-----------

#serie de tiempo de los errores
serieEt <- ts(residuals(mod1),freq=4,start=c(2000,1))
auto.arima(serieEt,ic="bic")
auto.arima(serieEt,ic="aic")
 
 

#-------------------------------------------ACF,PACF, LB, WD de los residuales-------------------------


acf(residuals(mod1),ci.type="ma",lag.max=24,lwd=2)
pacf(residuals(mod1),lag.max=24,lwd=2) 
BP.LB.test(residuals(mod1),maxlag=24,type="Ljung") 
pruebaDW1(mod1)

#------------------------------------------Seleccion de modelos de residuales---------------------------------------

eacf(residuals(mod1),ar.max=16 ,ma.max = 16)#Obteniendo las EACF
auto.arima(residuals(mod1),ic="bic")
auto.arima(residuals(mod1),ic="aic")
SelectModel(residuals(mod1),lag.max=36,Criterio="AIC",ARModel="AR")
SelectModel(residuals(mod1),lag.max=36,Criterio="BIC",ARModel="AR")

win.graph()
plot(armasubsets(residuals(mod1),nar=12,nma=12,y.name="AR",ar.method='ml'))
win.graph()
plot(armasubsets(residuals(mod1),nar=16,nma=16,y.name="AR",ar.method='ml'))
win.graph()
plot(armasubsets(residuals(mod1),nar=20,nma=20,y.name="AR",ar.method='ml'))

#-------------------------------------------Ajuste y pronostico--------------------


t=1:n
t2=t^2
t3=t^3
t4=t^4
indicadoras= seasonaldummy(yt)
X= cbind(t,t2,t3,t4,indicadoras)

#Valores de los predictores en los pron´osticos ex-post 
tnuevo=(n+1):(n+4) #son 4 pron´osticos
t2nuevo=tnuevo^2
t3nuevo=tnuevo^3
t4nuevo=tnuevo^4
trimestrenuevo=seasonaldummy(lnyt,h=4) #h=4 pues son cuatro pron´osticos despu´es
Xnuevo=cbind(t=tnuevo,t2=t2nuevo,t3=t3nuevo,t4nuevo,trimestre=trimestrenuevo) #esta matriz se usa en funci´on forecast
ytf=ts(Datos[tnuevo],frequency=4,start=c(2016,4)) #Los ´ultimos 4 valores de la serie, guardados para comparar


#------------------------------------AR(16)---------------------------------------------------------------------------------------------------------------------------
modelo1=Arima(lnyt,order=c(16,0,0),xreg=X,fixed=c(NA,0,0,0,NA,0,0,NA,
                0,NA,0,0,0,0,0,NA,rep(NA,8)),method="ML")
k1=24 #No. de par´ametros del modelo 2
df1=n-k1 #grados de libertad del modelo 2
coeftest(modelo1,df=df1)


ythat1=exp(modelo1$fitted)*exp(modelo1$sigma2/2) #valores ajustados en escala original
Res.origmodelo1=yt-ythat1 #seudo residuales
AICmodelo1=exp(crit.inf.resid(residuales= Res.origmodelo1,n.par=k1));AICmodelo1
BICmodelo1=exp(crit.inf.resid(residuales= Res.origmodelo1,n.par=k1,AIC="FALSE"));BICmodelo1


#An´alisis residuos de ajuste
resid.ajust1=residuals(modelo1) #residuos de ajuste
plot(resid.ajust1,main="Residuales vs. tiempo modelo1")
abline(h=0,lty=1)
abline(h=c(-2*sqrt(modelo1$sigma2),2*sqrt(modelo1$sigma2)),lty=2)
plot(as.numeric(modelo1$fitted),resid.ajust1,main="Residuales vs. ajustados modelo1")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo1$sigma2),2*sqrt(modelo1$sigma2)),lty=2)
acf(as.numeric(resid.ajust1),ci.type="ma",main="ACF modelo1",lag.max=24)
pacf(as.numeric(resid.ajust1),main="PACF modelo1",lag.max=24)
BP.LB.test(resid.ajust1,maxlag=24,type="Ljung")

#Evaluaci´on supuesto de normalidad
shapiro.test(resid.ajust1)
qqnorm(resid.ajust1,main="Gráfico Normal Residuales Modelo 1")
qqline(resid.ajust1,col=2,lwd=2)

#Para las predicciones se exponencia para llevar a escala original y se usa factor de correcci´on
predmod1=forecast(modelo1,xreg=Xnuevo,level=95) #predicciones en escala log
#Predicci´on en escala original
predic.modelo1=ts(exp(as.data.frame(predmod1))*exp(modelo1$sigma2/2),freq=4,start=c(2016,4))
predic.modelo1

#Calidad de pron´osticos
ytpron1=predic.modelo1[,1] #pron´osticos puntuales
accuracy(ytpron1,ytf)
Amplmodelo1=amplitud(LIP=predic.modelo1[,2],LSP=predic.modelo1[,3]); Amplmodelo1
Cobmodelo1=cobertura(real=ytf,LIP=predic.modelo1[,2],LSP=predic.modelo1[,3]); Cobmodelo1

#graficando la serie original, ajustada y sus pron´osticos, en escala original
plot(Datos)
lines(ythat1,col=2,lwd=2)
lines(ytpron1,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 1","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))


#----------------------------------ARMA(2,3)-------------------------------------------------------------------------------------------------------------------------
modelo2=Arima(lnyt,order=c(2,0,3),xreg=X,method="ML")
k2=13 #No. par´ametros modelo 5
df2=n-k2 #grados de libertad modelo 5
coeftest(modelo2,df=df2)

ythat2=exp(modelo2$fitted)*exp(modelo2$sigma2/2)
Res.origmodelo2=yt-ythat2 #seudo-residuos
AICmodelo2=exp(crit.inf.resid(residuales= Res.origmodelo2,n.par=k2)) ;AICmodelo2
BICmodelo2=exp(crit.inf.resid(residuales= Res.origmodelo2,n.par=k2,AIC="FALSE"));BICmodelo2
#An´alisis de residuales
resid.ajust2=residuals(modelo2) #residuos de ajuste
plot(resid.ajust2,main="Residuales vs. tiempo modelo2")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo5$sigma2),2*sqrt(modelo5$sigma2)),lty=2)
plot(as.numeric(modelo2$fitted),resid.ajust2,main="Residuales vs. ajustados modelo2")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo2$sigma2),2*sqrt(modelo2$sigma2)),lty=2)
acf(as.numeric(resid.ajust2),ci.type="ma",main="ACF modelo2",lag.max=24)
pacf(as.numeric(resid.ajust2),main="PACF modelo2",lag.max=24)
BP.LB.test(residuals(modelo2),maxlag=24,type="Ljung") 


#Evaluaci´on supuesto de normalidad
shapiro.test(resid.ajust2)
qqnorm(resid.ajust2,main="Gr´afico Normal Residuales Modelo 5")
qqline(resid.ajust2,col=2,lwd=2)

#Para las predicciones se exponencia para llevar a escala original y se aplica factor de correcci´on
predmod2=forecast(modelo2,xreg =Xnuevo,level=95) #predicciones en escala log
#Predicci´on en escala original
predic.modelo2=ts(exp(as.data.frame(predmod2))*exp(modelo2$sigma2/2),freq=4,start=c(2016,4))
predic.modelo2

#Calidad de pron´osticos
ytpron2=predic.modelo2[,1] #pron´osticos puntuales
accuracy(ytpron2,ytf)
Amplmodelo2=amplitud(LIP=predic.modelo2[,2],LSP=predic.modelo2[,3]); Amplmodelo2
Cobmodelo2=cobertura(real=ytf,LIP=predic.modelo2[,2],LSP=predic.modelo2[,3]); Cobmodelo2

#graficando la serie original, ajustada y sus pron´osticos, en escala original
plot(Datos)
lines(ythat2,col=2,lwd=2)
lines(ytpron2,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 2","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))

 

#----------------------------------ARMA(2,0)X(0,2)[4]-------------------------------------------------------------------------------------------------------------------
modelo3=Arima(lnyt,order=c(2,0,0),seasonal=list(order=c(0,0,2)),xreg=X,method="ML")
k3=12 #No. de par´ametros modelo 8
df3=n-k3 #grados de libertad modelo 8


ythat3=exp(modelo3$fitted)*exp(modelo3$sigma2/2)
Res.origmodelo3=yt-ythat3 #seudo-residuos
AICmodelo3=exp(crit.inf.resid(residuales= Res.origmodelo3,n.par=k3));AICmodelo3
BICmodelo3=exp(crit.inf.resid(residuales= Res.origmodelo3,n.par=k3,AIC="FALSE"));BICmodelo3

#An´alisis de residuales
resid.ajust3=residuals(modelo3) #residuos de ajuste
plot(resid.ajust3,main="Residuales vs. tiempo modelo3")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo3$sigma2),2*sqrt(modelo3$sigma2)),lty=2)
plot(as.numeric(modelo3$fitted),resid.ajust3,main="Residuales vs. ajustados modelo3")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo3$sigma2),2*sqrt(modelo3$sigma2)),lty=2)
acf(as.numeric(resid.ajust3),ci.type="ma",main="ACF modelo3",lag.max=24)
pacf(as.numeric(resid.ajust3),main="PACF modelo3",lag.max=24)
BP.LB.test(resid.ajust3,maxlag=24,type="Ljung")

#Evaluaci´on supuesto de normalidad
shapiro.test(resid.ajust3)
qqnorm(resid.ajust3,main="Gr´afico Normal Residuales Modelo 3")
qqline(resid.ajust3,col=2,lwd=2)

#Para las predicciones se exponencia para llevar a escala original y se aplica factor de correcci´on
predmod3=forecast(modelo3,xreg =Xnuevo,level=95) #predicciones en escala log
#Predicci´on destranformada a escala original
predic.modelo3=ts(exp(as.data.frame(predmod3))*exp(modelo3$sigma2/2),freq=4,start=c(2016,4))
predic.modelo3

#Calidad de pron´osticos
ytpron3=predic.modelo3[,1] #pron´osticos puntuales
accuracy(ytpron3,ytf)
Amplmodelo3=amplitud(LIP=predic.modelo3[,2],LSP=predic.modelo3[,3]); Amplmodelo3
Cobmodelo3=cobertura(real=ytf,LIP=predic.modelo3[,2],LSP=predic.modelo3[,3]); Cobmodelo3

#graficando la serie original, ajustada y sus pron´osticos, en escala original
plot(Datos)
lines(ythat3,col=2,lwd=2)
lines(ytpron3,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 3","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))

#----------------------------------AR(4)-------------------------------------------------------------------------------------------------------------------------
modelo4=Arima(lnyt,order=c(4,0,0),xreg=X,fixed=c(NA,0,0,NA,rep(NA,8)),method="ML")
k4=12 #No. de par´ametros del modelo 2
df1=n-k4 #grados de libertad del modelo 2

ythat4=exp(modelo4$fitted)*exp(modelo4$sigma2/2) #valores ajustados en escala original
Res.origmodelo4=yt-ythat4 #seudo residuales
AICmodelo4=exp(crit.inf.resid(residuales= Res.origmodelo4,n.par=k4));AICmodelo4
BICmodelo4=exp(crit.inf.resid(residuales= Res.origmodelo4,n.par=k4,AIC="FALSE"));BICmodelo4

#An´alisis residuos de ajuste
resid.ajust4=residuals(modelo4) #residuos de ajuste
plot(resid.ajust4,main="Residuales vs. tiempo modelo4")
abline(h=0,lty=1)
abline(h=c(-2*sqrt(modelo4$sigma2),2*sqrt(modelo4$sigma2)),lty=2)
plot(as.numeric(modelo4$fitted),resid.ajust4,main="Residuales vs. ajustados modelo4")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo4$sigma2),2*sqrt(modelo4$sigma2)),lty=2)
acf(as.numeric(resid.ajust4),ci.type="ma",main="ACF modelo4",lag.max=24)
pacf(as.numeric(resid.ajust4),main="PACF modelo4",lag.max=24)
BP.LB.test(resid.ajust4,maxlag=24,type="Ljung")

#Evaluaci´on supuesto de normalidad
shapiro.test(resid.ajust4)
qqnorm(resid.ajust4,main="Gráfico Normal Residuales Modelo 4")
qqline(resid.ajust4,col=2,lwd=2)

#Para las predicciones se exponencia para llevar a escala original y se usa factor de correcci´on
predmod4=forecast(modelo4,xreg=Xnuevo,level=95) #predicciones en escala log
#Predicci´on en escala original
predic.modelo4=ts(exp(as.data.frame(predmod4))*exp(modelo4$sigma2/2),freq=4,start=c(2016,4))
predic.modelo4

#Calidad de pron´osticos
ytpron4=predic.modelo4[,1] #pron´osticos puntuales
accuracy(ytpron4,ytf)
Amplmodelo4=amplitud(LIP=predic.modelo4[,2],LSP=predic.modelo4[,3]); Amplmodelo4
Cobmodelo4=cobertura(real=ytf,LIP=predic.modelo4[,2],LSP=predic.modelo4[,3]); Cobmodelo4

#graficando la serie original, ajustada y sus pron´osticos, en escala original
plot(Datos)
lines(ythat4,col=2,lwd=2)
lines(ytpron4,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 4","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))

#----------------------------------AR(2)-----------------------------------------------------------------------------------------------------------------------

modelo5=Arima(lnyt,order=c(2,0,0),xreg=X,method="ML")
k5=10 #No. de par´ametros del modelo              2
df5=n-k5 #grados de libertad del modelo 2

ythat5=exp(modelo5$fitted)*exp(modelo5$sigma2/2) #valores ajustados en escala original
Res.origmodelo5=yt-ythat5 #seudo residuales
AICmodelo5=exp(crit.inf.resid(residuales= Res.origmodelo5,n.par=k5));AICmodelo5
BICmodelo5=exp(crit.inf.resid(residuales= Res.origmodelo5,n.par=k5,AIC="FALSE"));BICmodelo5

#An´alisis residuos de ajuste
resid.ajust5=residuals(modelo5) #residuos de ajuste
plot(resid.ajust5,main="Residuales vs. tiempo modelo5")
abline(h=0,lty=1)
abline(h=c(-2*sqrt(modelo5$sigma2),2*sqrt(modelo5$sigma2)),lty=2)
plot(as.numeric(modelo5$fitted),resid.ajust5,main="Residuales vs. ajustados modelo5")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo5$sigma2),2*sqrt(modelo5$sigma2)),lty=2)
acf(as.numeric(resid.ajust5),ci.type="ma",main="ACF modelo5",lag.max=24)
pacf(as.numeric(resid.ajust5),main="PACF modelo5",lag.max=24)
BP.LB.test(resid.ajust5,maxlag=24,type="Ljung")

#Evaluaci´on supuesto de normalidad
shapiro.test(resid.ajust5)
qqnorm(resid.ajust5,main="Gráfico Normal Residuales Modelo 5")
qqline(resid.ajust5,col=2,lwd=2)

#Para las predicciones se exponencia para llevar a escala original y se usa factor de correcci´on
predmod5=forecast(modelo5,xreg=Xnuevo,level=95) #predicciones en escala log
#Predicci´on en escala original
predic.modelo5=ts(exp(as.data.frame(predmod5))*exp(modelo5$sigma2/2),freq=4,start=c(2016,4))
predic.modelo5

#Calidad de pron´osticos
ytpron5=predic.modelo5[,1] #pron´osticos puntuales
accuracy(ytpron5,ytf)
Amplmodelo5=amplitud(LIP=predic.modelo5[,2],LSP=predic.modelo5[,3]); Amplmodelo5
Cobmodelo5=cobertura(real=ytf,LIP=predic.modelo5[,2],LSP=predic.modelo5[,3]); Cobmodelo5

#graficando la serie original, ajustada y sus pron´osticos, en escala original
plot(Datos)
lines(ythat5,col=2,lwd=2)
lines(ytpron5,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 5","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))


#----------------------------------AR(10)----------------------------------------------------------------------------------------------------------------------


modelo6=Arima(lnyt,order=c(10,0,0),fixed=c(rep(0,4),NA,rep(0,4),NA,rep(NA,8)),xreg=X,method="ML")
k6=18 #No. de par´ametros del modelo 2
df6=n-k6 #grados de libertad del modelo 2

ythat6=exp(modelo6$fitted)*exp(modelo6$sigma2/2) #valores ajustados en escala original
Res.origmodelo6=yt-ythat6 #seudo residuales
AICmodelo6=exp(crit.inf.resid(residuales= Res.origmodelo6,n.par=k6));AICmodelo6
BICmodelo6=exp(crit.inf.resid(residuales= Res.origmodelo6,n.par=k6,AIC="FALSE"));BICmodelo6

#An´alisis residuos de ajuste
resid.ajust6=residuals(modelo6) #residuos de ajuste
plot(resid.ajust6,main="Residuales vs. tiempo modelo6")
abline(h=0,lty=1)
abline(h=c(-2*sqrt(modelo6$sigma2),2*sqrt(modelo6$sigma2)),lty=2)
plot(as.numeric(modelo6$fitted),resid.ajust6,main="Residuales vs. ajustados modelo6")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo6$sigma2),2*sqrt(modelo6$sigma2)),lty=2)
acf(as.numeric(resid.ajust6),ci.type="ma",main="ACF modelo6",lag.max=24)
pacf(as.numeric(resid.ajust6),main="PACF modelo6",lag.max=24)
BP.LB.test(resid.ajust6,maxlag=24,type="Ljung")

#Evaluaci´on supuesto de normalidad
shapiro.test(resid.ajust6)
qqnorm(resid.ajust6,main="Gráfico Normal Residuales Modelo 6")
qqline(resid.ajust6,col=2,lwd=2)

#Para las predicciones se exponencia para llevar a escala original y se usa factor de correcci´on
predmod6=forecast(modelo6,xreg=Xnuevo,level=95) #predicciones en escala log
#Predicci´on en escala original
predic.modelo6=ts(exp(as.data.frame(predmod6))*exp(modelo6$sigma2/2),freq=4,start=c(2016,4))
predic.modelo6

#Calidad de pron´osticos
ytpron6=predic.modelo6[,1] #pron´osticos puntuales
accuracy(ytpron6,ytf)
Amplmodelo6=amplitud(LIP=predic.modelo6[,2],LSP=predic.modelo6[,3]); Amplmodelo6
Cobmodelo6=cobertura(real=ytf,LIP=predic.modelo6[,2],LSP=predic.modelo6[,3]); Cobmodelo6

#graficando la serie original, ajustada y sus pron´osticos, en escala original
plot(Datos)
lines(ythat6,col=2,lwd=2)
lines(ytpron6,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 6","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))


#----------------------------------ARMA(5,7)--------------------------------------------------------------------------------------------------------------

modelo7=Arima(lnyt,order=c(5,0,7),fixed=c(NA,rep(0,3),NA,rep(0,6),NA,rep(NA,8)),xreg=X,method="ML")
k7=20 #No. de par´ametros modelo 6
df7=n-k7 #grados de libertad modelo 6

ythat7=exp(modelo7$fitted)*exp(modelo7$sigma2/2)
Res.origmodelo7=yt-ythat7 #seudo-residuos
AICmodelo7=exp(crit.inf.resid(residuales= Res.origmodelo7,n.par=k7)) ;AICmodelo7
BICmodelo7=exp(crit.inf.resid(residuales= Res.origmodelo7,n.par=k7,AIC="FALSE"));BICmodelo7

#An´alisis de residuales
resid.ajust7=residuals(modelo7) #residuos de ajuste
plot(resid.ajust7,main="Residuales vs. tiempo modelo7")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo7$sigma2),2*sqrt(modelo7$sigma2)),lty=2)
plot(as.numeric(modelo7$fitted),resid.ajust7,main="Residuales vs. ajustados modelo7")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo7$sigma2),2*sqrt(modelo7$sigma2)),lty=2)
acf(as.numeric(resid.ajust7),ci.type="ma",main="ACF modelo7",lag.max=24)
pacf(as.numeric(resid.ajust7),main="PACF modelo7",lag.max=24)
BP.LB.test(resid.ajust7,maxlag=24,type="Ljung")

#Evaluaci´on supuesto de normalidad
shapiro.test(resid.ajust7)
qqnorm(resid.ajust7,main="Gr´afico Normal Residuales Modelo 7")
qqline(resid.ajust7,col=2,lwd=2)

#Para las predicciones se exponencia para llevar a escala original y se aplica factor de correcci´on
predmod7=forecast(modelo7,xreg =Xnuevo,level=95) #predicciones en escala log
#Predicci´on en escala original
predic.modelo7=ts(exp(as.data.frame(predmod7))*exp(modelo7$sigma2/2),freq=4,start=c(2016,4))
predic.modelo7

#Calidad de pron´osticos
ytpron7=predic.modelo7[,1] #pron´osticos puntuales
accuracy(ytpron7,ytf)
Amplmodelo7=amplitud(LIP=predic.modelo7[,2],LSP=predic.modelo7[,3]); Amplmodelo7
Cobmodelo7=cobertura(real=ytf,LIP=predic.modelo7[,2],LSP=predic.modelo7[,3]); Cobmodelo7

#graficando la serie original, ajustada y sus pron´osticos, en escala original
plot(Datos)
lines(ythat7,col=2,lwd=2)
lines(ytpron7,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 7","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))

#------------------------------------------------ARIMA(2,1)X(4,3)[4]----------------------------------------------

modelo8=Arima(lnyt,order=c(2,0,1),seasonal=list(order=c(4,0,3)),
                                 xreg=X,method="CSS")
k8=18 #No. de par´ametros modelo 8
df8=n-k8 #grados de libertad modelo 8
coeftest(modelo8,df=df8)

ythat8=exp(modelo8$fitted)*exp(modelo8$sigma2/2)
Res.origmodelo8=yt-ythat8 #seudo-residuos
AICmodelo8=exp(crit.inf.resid(residuales= Res.origmodelo8,n.par=k8));AICmodelo8
BICmodelo8=exp(crit.inf.resid(residuales= Res.origmodelo8,n.par=k8,AIC="FALSE"));BICmodelo8

#An´alisis de residuales
resid.ajust8=residuals(modelo8) #residuos de ajuste
plot(resid.ajust8,main="Residuales vs. tiempo modelo3",xlim=c(2004,2017))
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo8$sigma2),2*sqrt(modelo8$sigma2)),lty=2)
plot(as.numeric(modelo8$fitted),resid.ajust8,xlim=c(11.1,12.1),main="Residuales vs. ajustados modelo3")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo8$sigma2),2*sqrt(modelo8$sigma2)),lty=2)
acf(as.numeric(resid.ajust8),ci.type="ma",main="ACF modelo3",lag.max=24)
pacf(as.numeric(resid.ajust8),main="PACF modelo3",lag.max=24)
BP.LB.test(resid.ajust8,maxlag=24,type="Ljung") 


#Evaluaci´on supuesto de normalidad
shapiro.test(resid.ajust8[19:n])
qqnorm(resid.ajust8[19:n],main="Gráfico Normal Residuales Modelo 3")
qqline(resid.ajust8[19:n],col=2,lwd=2)

#Para las predicciones se exponencia para llevar a escala original y se aplica factor de correcci´on
predmod3=forecast(modelo8,xreg =Xnuevo,level=95) #predicciones en escala log
#Predicci´on destranformada a escala original
predic.modelo8=ts(exp(as.data.frame(predmod3))*exp(modelo3$sigma2/2),freq=4,start=c(2016,4))
predic.modelo8

#Calidad de pron´osticos
ytpron8=predic.modelo8[,1] #pron´osticos puntuales
accuracy(ytpron8,ytf)
Amplmodelo8=amplitud(LIP=predic.modelo8[,2],LSP=predic.modelo8[,3]); Amplmodelo8
Cobmodelo8=cobertura(real=ytf,LIP=predic.modelo8[,2],LSP=predic.modelo8[,3]); Cobmodelo8

#graficando la serie original, ajustada y sus pron´osticos, en escala original
plot(Datos)
lines(ythat8,col=2,lwd=2)
lines(ytpron8,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 3","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))

#---------------------------ARMA(10,17) propuesto en asesoria-------------------------
modelo9=Arima(lnyt,order=c(10,0,17),xreg=X,fixed=c(0,NA,0,0,NA,rep(0,4),NA,0,0,0,NA,rep(0,7),NA,rep(0,4), NA,rep(NA,8)),method="ML")
k9=35 #No. de par´ametros del modelo 2
df9=n-k9 #grados de libertad del modelo 2

ythat9=exp(modelo9$fitted)*exp(modelo9$sigma2/2) #valores ajustados en escala original
Res.origmodelo9=yt-ythat9 #seudo residuales
AICmodelo9=exp(crit.inf.resid(residuales= Res.origmodelo9,n.par=k9));AICmodelo4
BICmodelo9=exp(crit.inf.resid(residuales= Res.origmodelo9,n.par=k9,AIC="FALSE"));BICmodelo4

#An´alisis residuos de ajuste
resid.ajust9=residuals(modelo9) #residuos de ajuste
plot(resid.ajust9,main="Residuales vs. tiempo modelo9")
abline(h=0,lty=1)
abline(h=c(-2*sqrt(modelo9$sigma2),2*sqrt(modelo9$sigma2)),lty=2)
plot(as.numeric(modelo9$fitted),resid.ajust9,main="Residuales vs. ajustados modelo9")
abline(h=0,lty=2)
abline(h=c(-2*sqrt(modelo9$sigma2),2*sqrt(modelo9$sigma2)),lty=2)
acf(as.numeric(resid.ajust9),ci.type="ma",main="ACF modelo9",lag.max=24)
pacf(as.numeric(resid.ajust9),main="PACF modelo9",lag.max=24)
BP.LB.test(resid.ajust9,maxlag=24,type="Ljung")

#Evaluaci´on supuesto de normalidad
shapiro.test(resid.ajust9)
qqnorm(resid.ajust9,main="Gráfico Normal Residuales Modelo 9")
qqline(resid.ajust9,col=2,lwd=2)

#Para las predicciones se exponencia para llevar a escala original y se usa factor de correcci´on
predmod9=forecast(modelo9,xreg=Xnuevo,level=95) #predicciones en escala log
#Predicci´on en escala original
predic.modelo9=ts(exp(as.data.frame(predmod4))*exp(modelo1$sigma2/2),freq=4,start=c(2016,4))
predic.modelo9

#Calidad de pron´osticos
ytpron9=predic.modelo9[,1] #pron´osticos puntuales
accuracy(ytpron9,ytf)
Amplmodelo9=amplitud(LIP=predic.modelo4[,2],LSP=predic.modelo4[,3]); Amplmodelo4
Cobmodelo9=cobertura(real=ytf,LIP=predic.modelo4[,2],LSP=predic.modelo4[,3]); Cobmodelo4

#graficando la serie original, ajustada y sus pron´osticos, en escala original
plot(Datos)
lines(ythat9,col=2,lwd=2)
lines(ytpron9,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 9","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))

#------------------------------------------estabilidad-----------------------------------------------------------

##Creando funci?on para estimaci?on recursiva del modelo
estim=function(n){
modelo1=lm(log(Datos)[1:n]~t[1:n]+I(t^2)[1:n]+I(t^3)[1:n]+I1[1:n]+I2[1:n]+I3[1:n]) 
resul1=cbind(coef(mod1),confint(mod1))
resul1 
}

#Definicion de tama?o de muestra para las regresiones recursivas
n=matrix(8:length(Datos),ncol=1) 
 

#Realizando las regresiones recursivas. En el arreglo se guardan matrices #cada una de 14 filas (se estiman p=7 par?ametros) y 3 columnas (estimaci?on, LIC, LSC)
b=array(apply(n,1,estim),dim=c(7,3,nrow(n))) 

betas0=b[1,,] #Extrayendo informaci?on de los interceptos
betas1=b[2,,] #Extrayendo informaci?on de los betas1
betas2=b[3,,] #Extrayendo informaci?on de los betas2
betas2=b[4,,] #Extrayendo informaci?on de los betas3
Ind_1=b[5,,] #Extrayendo informaci?on de los deltas1
Ind_2=b[6,,] #Extrayendo informaci?on de los deltas2
Ind_3=b[7,,] #Extrayendo informaci?on de los deltas3

betas0
betas1

rr1 <- recresid(mod1)
plot(rr1, type = "l",ylab="Residuales recursivos",xlab="t")
abline(h=0,col=2,lwd=2)
plot(efp(log(yt)~t+I(t^2)+I(t^3)+I(t^4)+I1+I2+I3, type = "Rec-CUSUM"),lwd=2,alpha = 0.05)
#Test CUSUM recursivo para cambio estructural
sctest(log(yt)~t+I(t^2)+I(t^3)+I(t^4)+I1+I2+I3,type="Rec-CUSUM")



#modelo exponencial polinomial


rr2 =recresid(mod2)
plot(rr2, type = "l",ylab="Residuales recursivos",xlab="t") 
abline(h=0,col=2,lwd=2)
bin <- coef(mod1)

b0=bin[1]
b1=bin[2]
b2=bin[3]
b3=bin[4]
d1=bin[5]
d2=bin[6]
d3=bin[7]

plot(efp(log(yt)~b0+b1*t+b2*I(t^2)+b3*I(t^3)+d1*I1+d2*I2+d3*I3),
,type = "Rec-CUSUM"),lwd=2)
sctest(log(yt)~(b0+b1*t+b2*I(t^2)+b3*I(t^3)+d1*I1+d2*I2+d3*I3),type="Rec-CUSUM")

#-------------------------------------pronostico---------------------------------------

plot(ytnuevo,main="Valores reales y pronósticos últimos Q=4 trimestres\nModelos 1,2,3",type="b",pch=1,lty=1,col=1,lwd=2,  
    ylab="vértice en miles de millones de pesos",ylim=c(min(ytnuevo,ytpron1,ytpron2),max(ytnuevo,ytpron1,ytpron2)),xaxt="n")
lines(ytpron1,col=2,pch=2,lty=2,type="b",lwd=2)
lines(ytpron2,col=3,pch=3,lty=3,type="b",lwd=2) 
lines(ytpron9,col=4,pch=4,lty=4,type="b",lwd=2)
legend("top",legend=c("Observado","Pron. modelo 1","Pron. modelo 2","Pron. modelo 3"),
bty="n",col=c(1,2,3,4),pch=1:5,lty=c(1,2,3,4),lwd=2) 
axis(1,at=time(ytnuevo),labels=c("Trimestre1","Trimestre2","Trimestre3","Trimestre4")) 
 




















 




