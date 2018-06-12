library(forecast)
library(TSA)
library(strucchange)
library(fANCOVA)
library(FitAR)
library(lmtest)
library(car)
library(pdR)

#----------------------------------------------------------------------

#Creando funci´on usuario crit.inf.resid() para calcular C*n(p)
crit.inf.resid=function(residuales,n.par,AIC="TRUE"){
if(AIC=="TRUE"){
#Calcula AIC
CI=log(mean(residualesˆ2))+2*n.par/length(residuales)
}
if(AIC=="FALSE"){
#Calcula BIC
CI=log(mean(residualesˆ2))+n.par*log(length(residuales))/length(residuales)
}
CI
}


#--------------------------------------------------------------------------------
Datos <- read.table(file.choose(),header=T,sep=";",skip=11,dec=",",colClasses=c(rep("NULL",11),"numeric",rep("NULL",6)))
Datos <- ts(Datos,freq=4,start=c(2000,1))
n <- length(Datos)-4
m <- round(n/4,0)
t <- 1:n
yt <- ts(Datos[t],freq=4,start=c(2000,1)) #valores de la serie en muestra de ajuste
lnyt <- log(yt)
#------------------------------Grafica serie-------------------------------------------------------------------------------

g1 <- plot(Datos,main="Serie original",ylim=c(min(Datos),max(Datos)))
g2 <- plot(lnyt,main="Logarítmo  serie",ylim=c(min(trans),max(trans)))

difd1 <- diff(lnyt)
plot(difd1,main=expression(paste(nabla,sep="",log(Y[t]))))
abline(h=mean(difd1))


difs4 <- diff(lnyt,lag=4)
plot(difs4,main=expression(paste(nabla[4],sep="",log(Y[t]))))
abline(h=mean(difs4))


difdD4 <- diff(diff(lnyt),lag=4)) 
plot(difdD4,main=expression(paste(nabla[4],sep="",nabla,sep="",log(Y[t]))))
abline(h=mean(difdD4))
abline(h=0,col=2)


acf(as.numeric(lnyt),lag.max=24,ci.type="ma",ci.col=2,col=4)
acf(as.numeric(difd1),lag.max=24,ci.type="ma",ci.col="red",main=expression(paste("ACF",sep=" ","de",sep=" ",nabla,sep="",log(Y[t]))))
acf(as.numeric(difs4),lag.max=24,ci.type="ma",ci.col="red")
acf(as.numeric(difdD4),lag.max=24,ci.type="ma",ci.col="red",lwd=3,main=expression(paste("ACF",sep=" ","de",sep=" ",nabla[4],nabla,sep="",log(Y[t]))))
abline(v=seq(4,24,4),col=4,lty=2)
win.graph()
pacf(as.numeric(difdD4),lag.max=24,ci.col="red",lwd=3,main=expression(paste("PACF",sep=" ","de",sep=" ",nabla[4],nabla,sep="",log(Y[t]))),ylim=c(-0.45,0.45))
abline(v=seq(4,24,4),col=4,lty=2)

#-------------------------------------------------------------------------
auto.arima(lnyt,ic="aic",seasonal.test="ocsb")
auto.arima(lnyt,ic="aic",seasonal.test="ch")
auto.arima(lnyt,ic="bic",seasonal.test="ocsb")
auto.arima(lnyt,ic="bic",seasonal.test="ch")




#--------------------------------------------------------------------------
plot(armasubsets(difdD4,nar=12,nma=12))

plot(armasubsets(difdD4,nar=16,nma=16))

plot(armasubsets(difdD4,nar=20,nma=20))





#-------------------------------------------------------------------------
HEGY.test(wts=difdD4,itsd=c(0,0,c(0)),selectlags=list(mode="aic", Pmax=4))$stats

#-------------------------------------------Ajuste y pronostico--------------------




#Valores de los predictores en los pron?osticos ex-post 
tnuevo=(n+1):(n+4) #son 4 pron?osticos

ytf=ts(Datos[tnuevo],frequency=4,start=c(2016,4)) #Los ?ultimos 4 valores de la serie, guardados para comparar




#--------------------------------------------------------------------------------------------
#MODELO 1 ARIMA(0,1,0)(0,1,1)[4]

modelo.1=Arima(lnyt,order=c(0,1,0),seasonal=list(order=c(0,1,1)),method="ML")
coeftest(modelo.1) #Tabla de par?metros estimados; valores P bajo la N(0,1)




#ACF Y PACF 
acf(as.numeric(residuals(modelo.1)),ci.type="ma",lag.max=24,lwd=2,main="ACF modelo 1")
pacf(as.numeric(residuals(modelo.1)),lag.max=24,lwd=2,main="PACF modelo 1")

#TEST LJUNG BOX 
resid.ajust1=residuals(modelo.1)
BP.LB.test(resid.ajust1,maxlag=24,type="Ljung")

#Test de Normalidad Shapiro y gr?fica
shapiro.test(residuals(modelo.1))
qqnorm(residuals(modelo.1),main="Grafíco normal residuales modelo 1");qqline(residuals(modelo.1),col=2)

#Grafica de Residuales vs Tiempo
plot((residuals(modelo.1)),main="Residuales vs Tiempo modelo 1");abline(h=0,col=2)
abline(h=c(-2*sqrt(modelo.1$sigma2),mean(residuals(modelo.1)),2*sqrt(modelo.1$sigma2)),lty=2,col=1)

#Grafica de Residuales vs Ajustados 
plot(as.numeric(modelo.1$fitted),residuals(modelo.1),main="Residuales vs Ajustados modelo1")
abline(h=c(-2*sqrt(modelo.1$sigma2),0,2*sqrt(modelo.1$sigma2)),lty=2,col=2)



k1=1#No. de par?ametros del modelo 2
df1=n-1 #grados de libertad del modelo 2


ythat1=exp(modelo.1$fitted)*exp(modelo.1$sigma2/2) #valores ajustados en escala original
Res.origmodelo1=yt-ythat1 #seudo residuales
AICmodelo1=exp(crit.inf.resid(residuales= Res.origmodelo1,n.par=k1));AICmodelo1
BICmodelo1=exp(crit.inf.resid(residuales= Res.origmodelo1,n.par=k1,AIC="FALSE"));BICmodelo1





pronm1=exp(as.data.frame(forecast(modelo.1,h=4,level=95)))*exp(modelo.1$sigma2/2) #pronósticos en escala original
pronm1=ts(pronm1,freq=4,start=c(2016,4))
pronm1
accuracy(pronm1[,1],ytf) 
Cobmodelo1=cobertura(real=ytf,LIP=pronm1[,2],LSP=pronm1[,3]); Cobmodelo1
Amplmodelo1=amplitud(LIP=pronm1[,2],LSP=pronm1[,3]); Amplmodelo1

#graficando la serie original, ajustada y sus pron?osticos, en escala original
plot(Datos)
lines(ythat1,col=2,lwd=2)
lines(pronm1[,1],col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 1"),lty=1,lwd=c(1,2),col=c(1,2))



#----------------------------------------------------------------
#//////////////////////////////////////////////////////////

#MODELO 2 ARIMA(0,1,0)(2,1,1)[4]

modelo.2=Arima(log(yt),order=c(0,1,0),seasonal=list(order=c(2,1,1)),method="ML")

coeftest(modelo.2) #Tabla de par?ametros estimados; valores P bajo la N(0,1)

#ACF Y PACF 
acf(as.numeric(residuals(modelo.2)),ci.type="ma",lag.max=24,lwd=2)

pacf(as.numeric(residuals(modelo.2)),lag.max=24,lwd=2)

#TEST LJUNG BOX 
resid.ajust2=residuals(modelo.2)
BP.LB.test(Res.origmodelo2,maxlag=24,type="Ljung")


#Test de Normalidad Shapiro y gr?fica
shapiro.test(residuals(modelo.2))
qqnorm(residuals(modelo.2));qqline(residuals(modelo.2),col=2)

#Grafica de Residuales
plot((residuals(modelo.2)));abline(h=0,col=1)
abline(h=c(-2*sqrt(modelo.2$sigma2),mean(residuals(modelo.2)),2*sqrt(modelo.2$sigma2)),lty=2,col=2)

#Grafica de reidauloes v
plot(as.numeric(modelo.2$fitted),residuals(modelo.2))
abline(h=c(-2*sqrt(modelo.2$sigma2),0,2*sqrt(modelo.2$sigma2)),lty=2,col=2)



k2=3#No. de par?ametros del modelo 2
df2=n-3 #grados de libertad del modelo 2


ythat2=exp(modelo.2$fitted)*exp(modelo.2$sigma2/2) #valores ajustados en escala original
Res.origmodelo2=yt-ythat2 #seudo residuales
AICmodelo2=exp(crit.inf.resid(residuales= Res.origmodelo2,n.par=k2));AICmodelo2
BICmodelo2=exp(crit.inf.resid(residuales= Res.origmodelo2,n.par=k2,AIC="FALSE"));BICmodelo2




pronm2=exp(as.data.frame(forecast(modelo.2,h=4,level=95)))*exp(modelo.2$sigma2/2) #pronósticos en escala original
pronm2=ts(pronm2,freq=4,start=c(2016,4))
pronm2

accuracy(pronm2[,1],ytf) 

Cobmodelo2=cobertura(real=ytf,LIP=pronm2[,2],LSP=pronm2[,3]); Cobmodelo2
Amplmodelo2=amplitud(LIP=pronm2[,2],LSP=pronm2[,3]); Amplmodelo2

#graficando la serie original, ajustada y sus pron?osticos, en escala original
plot(Datos)
lines(ythat2,col=2,lwd=2)
lines(pronm2[,1],col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 2","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))



//////////////////////////////////////
#-----------------------------------------------------------------------
#MODELO 3 ARIMA(3,1,0)(2,1,3)


modelo.3=Arima(log(yt),order=c(3,1,0),seasonal=list(order=c(2,1,3)),method="ML",
fixed=c(NA,0,NA,NA,NA,0,NA,NA))
coeftest(modelo.3)

#ACF Y PACF 
acf(as.numeric(residuals(modelo.3)),ci.type="ma",lag.max=24,lwd=2)

pacf(as.numeric(residuals(modelo.3)),lag.max=24,lwd=2)

#TEST LJUNG BOX 
resid.ajust3=residuals(modelo.3)
BP.LB.test(Res.origmodelo3,maxlag=24,type="Ljung")


#Test de Normalidad Shapiro y gr?fica
shapiro.test(residuals(modelo.3))
qqnorm(residuals(modelo.3));qqline(residuals(modelo.3),col=2)

#Grafica de Residuales
plot((residuals(modelo.3)));abline(h=0,col=1)
abline(h=c(-2*sqrt(modelo.3$sigma2),mean(residuals(modelo.3)),2*sqrt(modelo.3$sigma2)),lty=2,col=2)

#Grafica de reidauloes v
plot(as.numeric(modelo.3$fitted),residuals(modelo.3))
abline(h=c(-2*sqrt(modelo.3$sigma2),0,2*sqrt(modelo.3$sigma2)),lty=2,col=2)


k3=6#No. de par?ametros del modelo 2
df3=n-6 #grados de libertad del modelo 2


ythat3=exp(modelo.3$fitted)*exp(modelo.3$sigma2/2) #valores ajustados en escala original
Res.origmodelo3=yt-ythat3 #seudo residuales
AICmodelo3=exp(crit.inf.resid(residuales= Res.origmodelo3,n.par=k3));AICmodelo3
BICmodelo3=exp(crit.inf.resid(residuales= Res.origmodelo3,n.par=k3,AIC="FALSE"));BICmodelo3




pronm3=exp(as.data.frame(forecast(modelo.3,h=4,level=95)))*exp(modelo.3$sigma2/2) #pronósticos en escala original
pronm3=ts(pronm3,freq=4,start=c(2016,4))
pronm3

accuracy(pronm3[,1],ytf) 
Amplmodelo3=amplitud(LIP=pronm3[,2],LSP=pronm3[,3]); Amplmodelo3

Cobmodelo3=cobertura(real=ytf,LIP=pronm3[,2],LSP=pronm3[,3]); Cobmodelo3

#graficando la serie original, ajustada y sus pron?osticos, en escala original
plot(Datos)
lines(ythat3,col=2,lwd=2)
lines(pronm3[,1],col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 3","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))



#///////////////////////////////////
#------------------------------------------------------------
#MODELO 4 ARIMA(3,1,3)(3,1,4)[4]


modelo.4=Arima(log(yt),order=c(3,1,3),seasonal=list(order=c(3,1,4)),method="ML",
fixed=c(0,NA,NA,NA,0,NA,0,NA,NA,NA,0,0,NA))

coeftest(modelo.4)


#ACF Y PACF 
acf(as.numeric(residuals(modelo.4)),ci.type="ma",lag.max=24,lwd=2,main="ACF modelo2")

pacf(as.numeric(residuals(modelo.4)),lag.max=24,lwd=2,main="PACF modelo 2")

#TEST LJUNG BOX 
resid.ajust4=residuals(modelo.4)
BP.LB.test(Res.origmodelo4,maxlag=24,type="Ljung")

#Test de Normalidad Shapiro y gr?fica
shapiro.test(residuals(modelo.4))
qqnorm(residuals(modelo.4),main="Grafíco normal residuales modelo 2");qqline(residuals(modelo.4),col=2)



#Grafica de Residuales
plot((residuals(modelo.4)),main="Residuales vs Tiempo modelo 2",ylab="Residuales");abline(h=0,col=1)
abline(h=c(-2*sqrt(modelo.4$sigma2),mean(residuals(modelo.4)),2*sqrt(modelo.4$sigma2),lty=2,col=2)

#Grafica de reidauloes v
plot(as.numeric(modelo.4$fitted),residuals(modelo.4),main="Residuales vs Ajustados modelo2",ylab="Residuales")
abline(h=c(-2*sqrt(modelo.4$sigma2),0,2*sqrt(modelo.4$sigma2)),lty=2,col=2)



k4=8 #No. de par?ametros del modelo 2
df4=n-8 #grados de libertad del modelo 2


ythat4=exp(modelo.4$fitted)*exp(modelo.4$sigma2/2) #valores ajustados en escala original
Res.origmodelo4=yt-ythat4 #seudo residuales
AICmodelo4=exp(crit.inf.resid(residuales= Res.origmodelo4,n.par=k4));AICmodelo4
BICmodelo4=exp(crit.inf.resid(residuales= Res.origmodelo4,n.par=k4,AIC="FALSE"));BICmodelo4




pronm4=exp(as.data.frame(forecast(modelo.4,h=4,level=95)))*exp(modelo.4$sigma2/2) #pronósticos en escala original
pronm4=ts(pronm4,freq=4,start=c(2016,4))
pronm4

accuracy(pronm4[,1],ytf) 

Amplmodelo4=amplitud(LIP=pronm4[,2],LSP=pronm4[,3]); Amplmodelo4
Cobmodelo4=cobertura(real=ytf,LIP=pronm4[,2],LSP=pronm4[,3]); Cobmodelo4

#graficando la serie original, ajustada y sus pron?osticos, en escala original
plot(Datos,main="Serie real y ajustada\nAjuste SARIMA(3,1,3)(3,2,4)[4]")

plot(Datos)
lines(ythat4,col=2,lwd=2)
lines(pronm4[,1],col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado modelo 2"),lty=1,lwd=c(1,2),col=c(1,2))

 
#--------------------------------------------------------------

holt.winter <- HoltWinters(yt,gamma=0.9999,seasonal="multiplicative")

exp(crit.inf.resid(residuales=residuals(holt.winter),n.par=5))
exp(crit.inf.resid(residuales=residuals(holt.winter),n.par=5,AIC="FALSE"))

pronholt=predict(holt.winter,n.ahead=4,prediction=T,level = 0.95)
pronholt
ytpronholt=pronholt[,1]
accuracy(ytpronholt,ytnuevo)
AmplmodeloHol=amplitud(LIP=pronholt[,2],LSP=pronholt[,3]); AmplmodeloHol
CobmodeloHol=cobertura(real=ytf,LIP=pronholt[,2],LSP=pronholt[,3]); CobmodeloHol

plot(Datos,main="Serie real y ajustada\nAjuste por Suavizamiento Holt-Winters")
lines(fitted(holt.winter)[,1],col=2) 
lines(ytpronholt,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))


leg


#--------------------------------------------------------------------

t=1:n
t2=t^2
t3=t^3
t4=t^4
indicadoras= seasonaldummy(yt)
X= cbind(t,t2,t3,t4,indicadoras)

#Valores de los predictores en los pron?osticos ex-post 
tnuevo=(n+1):(n+4) #son 4 pron?osticos
t2nuevo=tnuevo^2
t3nuevo=tnuevo^3
t4nuevo=tnuevo^4
trimestrenuevo=seasonaldummy(lnyt,h=4) #h=4 pues son cuatro pron?osticos despu?es
Xnuevo=cbind(t=tnuevo,t2=t2nuevo,t3=t3nuevo,t4nuevo,trimestre=trimestrenuevo) #esta matriz se usa en funci?on forecast
ytf=ts(Datos[tnuevo],frequency=4,start=c(2016,4)) #Los ?ultimos 4 valores de la serie, guardados para comparar



#--------------------------------------------------------------



modelo2 <- Arima(log(yt),order=c(4,0,0),xreg=X,fixed=c(NA,rep(0,2),NA,rep(NA,8)),method="ML")


#Para las predicciones se exponencia para llevar a escala original y se aplica factor de correcci?on
predmod2=forecast(modelo2,xreg =Xnuevo,level=95) #predicciones en escala log
#Predicci?on en escala original
predic.modelo2=ts(exp(as.data.frame(predmod2))*exp(modelo2$sigma2/2),freq=4,start=c(2016,4))
predic.modelo2

#Calidad de pron?osticos
ytpron2=predic.modelo2[,1] #pron?osticos puntuales
accuracy(ytpron2,ytf)
Amplmodelo2=amplitud(LIP=predic.modelo2[,2],LSP=predic.modelo2[,3]); Amplmodelo2
Cobmodelo2=cobertura(real=ytf,LIP=predic.modelo2[,2],LSP=predic.modelo2[,3]); Cobmodelo2

ythat2=exp(modelo2$fitted)*exp(modelo2$sigma2/2)

plot(Datos,main="Serie real y ajustada\n modelo con errores AR(4)")
lines(ythat2,col=2,lwd=2)
lines(ytpron2,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))



#-----------------------------------------------------------------------------

t <- 1:n
yt <- ts(Datos[t],freq=4,start=c(2000,1)) #valores de la serie en muestra de ajuste
lnyt <- log(yt)
trimestre <- seasonaldummy(yt) #Matriz con las 3 primeras variables Indicadoras 
#Separando una a una las 3 variables indicadoras
I1 <- trimestre[,1]
I2 <- trimestre[,2]
I3 <- trimestre[,3]
#Valores ajustados

#-----------------------------------------------------------------------------

#Modelo log-polinomial estacional
mod1 <- lm(lnyt~t+I(t^2)+I(t^3)+I(t^4)+I1+I2+I3)


tnuevo <- (n+1):length(Datos) #?ndice de tiempo en los pron?sticos
ytnuevo=ts(Datos[tnuevo],freq=4,start=c(2016,4)) #valores reales de la serie en #muestra de validaci?on cruzada
trimestre.nuevo=seasonaldummy(yt,h=4) #matriz con las 11 
I1n <- trimestre.nuevo[,1]
I2n <- trimestre.nuevo[,2]
I3n <- trimestre.nuevo[,3]


 

#Pronostico log polinomial
pronmod1=exp(predict(mod1,newdata=data.frame(t=tnuevo,I1=I1n,I2=I2n,I3=I3n)))*exp(summary(mod1)$sigma^2/2)
pronmod1=ts(pronmod1,freq=4,start=c(2016,4))
pronmod1
accuracy(pronmod1,ytnuevo)
yhatmod1=ts(exp(fitted(mod1))*exp(summary(mod1)$sigma^2/2),freq=4,start=c(2000,1))

plot(Datos,main="Serie real y ajustada\n modelo de grado cuatro con indicadoras")
lines(yhatmod1,col=2,lwd=2)
lines(pronmod1,col="blue",lwd=2)
legend("topleft",legend=c("Original","ajustado","pronosticado"),lty=1,lwd=c(1,2,2),col=c(1,2,4))



#-----------------------------------------------------------------------------


plot(ytnuevo,main="Valores reales y pronósticos últimos\n Q=4 trimestres",type="b",pch=1,lty=1,col=1,lwd=2,  
    ylab="vértice en miles de millones de pesos",ylim=c(min(ytnuevo,pronm4[,1],pronm1[,1]),max(ytnuevo,pronm4[,1],pronm1[,1])),xaxt="n")
lines(pronm4[,1],col=2,pch=2,lty=2,type="b",lwd=2)
lines(pronm1[,1],col=3,pch=3,lty=3,type="b",lwd=2) 
legend("topright",legend=c("Observado","Pron. modelo SARIMA(3,1,3)(3,2,4)[4]","Pron. modelo SARIMA(0,1,0)(0,1,1)[4]"),
bty="n",col=c(1,2,3),pch=1:5,lty=c(1,2,3),lwd=2) 
axis(1,at=time(ytnuevo),labels=c("Trimestre1","Trimestre2","Trimestre3","Trimestre4")) 

#---------------------------------------------------------------------------------------------------------------------
win.graph()
plot(ytnuevo,main="Valores reales y pronósticos últimos\n Q=4 trimestres",pch=1,lty=1,col=1,lwd=2,  
    ylab="vértice en miles de millones de pesos",ylim=c(min(ytnuevo,pronmod1,predic.modelo2,pronm4[,1],ytpronholt),max(ytnuevo,pronmod1,predic.modelo2,pronm4,ytpronholt)),xaxt="n")
lines(pronmod1,col=2,pch=2,lty=2,type="b",lwd=2)
lines(ytpron2,col=3,pch=3,lty=3,type="b",lwd=2) 
lines(pronm4[,1],col=4,pch=4,lty=4,type="b",lwd=2) 
lines(ytpronholt,col=5,pch=5,lty=5,type="b",lwd=2) 
legend("topright",legend=c("Observado","Pron. modelo de grado 4 con indicadoras","Pron. modelos con errores AR(4)","Pron. SARIMA(3,1,3)(3,2,4)[4]","Pron. Holt Winters"),
bty="n",col=c(1,2,3,4,5),pch=1:5,lty=c(1,2,3,4,5),lwd=2) 
axis(1,at=time(ytnuevo),labels=c("Trimestre1","Trimestre2","Trimestre3","Trimestre4")) 
 





 



