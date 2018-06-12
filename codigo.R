library(forecast)
library(TSA)
library(strucchange)
library(fANCOVA)
Datos <- read.table(file.choose(),header=T,sep=";",skip=11,dec=",",colClasses=c(rep("NULL",11),"numeric",rep("NULL",6)))
Datos <- ts(Datos,freq=4,start=c(2000,1))
length(Datos)
######

#Creando funci?on para extraer correctamente los valores deltai
factoresdeltai=function(descom,s,estacionini){ if(estacionini==1){ deltasi=descom$figure }
  if(estacionini!=1){ j=estacionini;deltasi=c(descom$figure[(s-j+2):s],descom$figure[1:(s-j+1)]) } 
  deltasi
}

#Creando funci?on usuario crit.inf.resid() para calcular C_n?*(p)
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


#Punto 3.2

g1 <- plot(Datos,main="Serie original",ylim=c(min(Datos),max(Datos)))
trans <- log(Datos)
g2 <- plot(trans,main="Logar?tmo  serie",ylim=c(min(trans),max(trans)),xlab="t (Figura 2.2) ")

Tendencia <- decompose(log(Datos),type = "additive")$trend
plot(Tendencia,main="Tendencia de log(Yt)",xlab="t  ")
boxplot(trans~cycle(trans),main="Boxplot",xlab="t  ")
x=diff(trans)

periodogram(x,main="Periodograma",xlab = "F (Figura2.6)");abline(v=c(0.25,0.5),col=2,lty=2)

#Punto 3.4


n <- length(Datos)-4
t <- 1:n
yt <- ts(Datos[t],freq=4,start=c(2000,1)) #valores de la serie en muestra de ajuste
trimestre <- seasonaldummy(yt) #Matriz con las 11 primeras variables Indicadoras 
#Separando una a una las 11 variables indicadoras
I1 <- trimestre[,1]
I2 <- trimestre[,2]
I3 <- trimestre[,3]
#Ajuste de los modelos
#Modelo log-polinomial estacional
mod1 <- lm(log(yt)~t+I(t^2)+I(t^3)+I1+I2+I3)

#Guardando coeficientes estimados 
bin <- coef(mod1)

#Modelo exponencial cubico estacional
mod2 <- nls(yt~exp(b0+b1*t+b2*I(t^2)+b3*I(t^3)+d1*I1+d2*I2+d3*I3),
start =list(b0=bin[1],b1=bin[2],b2=bin[3],b3=bin[4],d1=bin[5],d2=bin[6],d3=bin[7]))

yhatmod1=ts(exp(fitted(mod1))*exp(summary(mod1)$sigma^2/2),freq=4,start=c(2000,1))#ajustado modelo 2

#C?lculo valores ajustados del modelo 2
yhatmod2=ts(fitted(mod2),freq=4,start=c(2000,1))
summary(mod1)#Resumen modelo 2

summary(mod2)#Resumen modelo 1


 anova(mod2)
#Grafica modelo 1
plot(Datos,main="Serie real y ajustada modelo1")
lines(yhatmod1,col=2,lwd=2)
legend("topleft",legend=c("Original","Ajustado modelo 1"),lty=1,col=c(1,2))

#Grafica modelo 2
plot(Datos,main="Serie real y ajustada modelo 2")
lines(yhatmod2,col=4,lwd=2)
legend("topleft",legend=c("Original"," Ajustado modelo 2"),lty=1,col=c(1,4))

#Suavizamiento exponencial Holt
holt.winter <- HoltWinters(yt,seasonal = "multiplicative",optim.start = c(alpha = 0.3, beta = 0.01, gamma = 0.1) )
holt.winter
#grafica holt
plot(Datos,main="Serie real y ajustada\nAjuste por Suavizamiento Holt-Winters")
lines(fitted(holt.winter)[,1],col=2) 
legend("topleft",legend=c("Original","Holt-Winters"),lty=1,col=c(1,2))

#Loess
yt.adj=yt/St
descom=decompose(yt,type="multiplicative")
St=descom$seasonal 
ajusteLoess=loess.as(t,yt.adj,degree=1,criterion="gcv",family="gaussian",plot=F) 
Tt=ts(fitted(ajusteLoess),frequency=4,start=c(2000,1)) #Tendencia ajustada LOESS lineal
ythat=Tt*St 
summary(ajusteLoess)

#Grafica loess

plot(Datos,main="Serie real y su ajuste\npor descomposici?n & LOESS lineal") 
lines(ythat,col=4)
legend("topleft",legend=c("Original","Ajustada"),col=c(1,4),lty=1) 

#deltas i loess
deltas_i=factoresdeltai(descom=descom,s=12,estacionini=1) #Obteniendo valor de los s 
data.frame(deltas_i)



#Calculo AIC y BIC aproximados version exp(C*n(p))
et1 <- yt-yhatmod1
et<- yt-ythat
p <- 5
p1

p1=round(ajusteLoess$enp)+3 

AIC2=exp(crit.inf.resid(residuales=et1,n.par=7))
BIC2=exp(crit.inf.resid(residuales=et1,n.par=7,AIC="FALSE"))

AIC1=exp(crit.inf.resid(residuales=residuals(mod2),n.par=7))
BIC1=exp(crit.inf.resid(residuales=residuals(mod2),n.par=7,AIC="FALSE"))

AIC4=exp(crit.inf.resid(residuales=et,n.par=p1))
BIC4=exp(crit.inf.resid(residuales=et,n.par=p1,AIC="FALSE"))

AIC3=exp(crit.inf.resid(residuales=residuals(holt.winter),n.par=p))
BIC3=exp(crit.inf.resid(residuales=residuals(holt.winter),n.par=p,AIC="FALSE"))

tablacriter=cbind(Log=c(AIC2,BIC2),expo=c(AIC1,BIC1),hol.win=c(AIC3,BIC3),loes=c(AIC4,BIC4))
rownames(tablacriter)=c("AIC","BIC")
tablacriter 


#3.5


#residuales modelo 1

plot.ts(residuals(mod1),ylim=c(min(residuals(mod1),-2*summary(mod1)$sigma),
max(residuals(mod1),2*summary(mod1)$sigma)),
main="Residuos vs. t\nModelo Log c?bico estacional")
abline(h=0,col=2)
abline(h=c(-2*summary(mod1)$sigma,2*summary(mod1)$sigma),col=2)
win.graph()
plot(fitted(mod1),residuals(mod1),ylim=c(min(residuals(mod1),-2*summary(mod1)$sigma),
max(residuals(mod1),2*summary(mod1)$sigma)),
main="Residuos vs. ajustados\nModelo log c?bico estacional")
abline(h=0,col=2)
abline(h=c(-2*summary(mod1)$sigma,2*summary(mod1)$sigma),col=2)

#residuales modelo 2


plot.ts(residuals(mod2),ylim=c(min(residuals(mod2),-2*summary(mod2)$sigma),
max(residuals(mod2),2*summary(mod2)$sigma)),
main="Residuos vs. t\nModelo exponencial c?bico estacional")
abline(h=0,col=2)
abline(h=c(-2*summary(mod2)$sigma,2*summary(mod2)$sigma),col=2)
win.graph()
plot(fitted(mod2),residuals(mod2),ylim=c(min(residuals(mod2),-2*summary(mod2)$sigma),
max(residuals(mod2),2*summary(mod2)$sigma)),
main="Residuos vs. ajustados\nModelo exponencial c?bico estacional")
abline(h=0,col=2)
abline(h=c(-2*summary(mod2)$sigma,2*summary(mod2)$sigma),col=2)

#residuales loess

df=n-7
MSE=sum(et^2)/df #MSE aproximado del ajuste total del modelo 1
plot(et,ylim=c(min(-2*sqrt(MSE),et),max(2*sqrt(MSE),et)),main="Residuos vs. t\nAjuste por descomposici?n & LOESS
lineal")
abline(h=0,col=2)
abline(h=c(-2*sqrt(MSE),2*sqrt(MSE)),col=2)
win.graph()
plot(as.numeric(ythat),et,ylim=c(min(-2*sqrt(MSE),et),max(2*sqrt(MSE),et)),main="Residuos vs. ajustados\nAjuste
por descomposici?n & LOESS lineal",xlab="ythat")
abline(h=0,col=2)
abline(h=c(-2*sqrt(MSE),2*sqrt(MSE)),col=2)

#residuales holt

df3=n-2*4-((4-1)+2)
MSE3=holt.winter$SSE/df3 #MSE aproximado del ajuste total del Suavizamiento
MSE3 
plot(residuals(holt.winter),ylim=c(min(-2*sqrt(MSE3),residuals(holt.winter)),max(2*sqrt(MSE3),residuals(holt.winter))),
 main="Residuos vs. t\nAjuste por suavizamiento exponencial Holt-Winters")
abline(h=0,col=2)
abline(h=c(-2*sqrt(MSE3),2*sqrt(MSE3)),col=2)
win.graph()
plot(as.numeric(fitted(holt.winter)[,1]),residuals(holt.winter),
 ylim=c(min(-2*sqrt(MSE3),residuals(holt.winter)),max(2*sqrt(MSE3),residuals(holt.winter))),
 main="Residuos vs. ajustados\nAjuste por suavizamiento exponencial Holt-Winters",xlab="fitted")
abline(h=0,col=2)
abline(h=c(-2*sqrt(MSE3),2*sqrt(MSE3)),col=2)

#3.6


#Valores de las variables en la validaci?n cruzada
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

#pronostico exponencial cubico estacionario
pronmod2=predict(mod2,newdata=data.frame(t=tnuevo,I1=I1n,I2=I2n,I3=I3n),interval="prediction",level=0.95) 
pronmod2=ts(pronmod2,freq=4,start=c(2016,4))
accuracy(pronmod2,ytnuevo)


#Pronostico holt
pronholt=predict(holt.winter,n.ahead=4,prediction=T,level = 0.95)
pronholt
ytpronholt=pronholt[,1]
accuracy(ytpronholt,ytnuevo)
dim(pronholt)


#Pronostico loess
alfa.optim=ajusteLoess$pars$span #guardando el valor ?optimo del par?ametro alfa
#Pron?osticos para la componente estacional
i=c(4,1,2,3) #identificando la estaci?on correspondiente a los m=4 per?iodos de pron?ostico
Stnuevo=deltas_i[i] #Asignando el valor de St a los per?iodos a pronosticar
Stnuevo=ts(Stnuevo,frequency=4,start=c(2016,4)) #convirtiendo en serie de tiempo al pron?ostico de St
#pronostico tendencia
Ttnuevo=predict(loess(yt.adj~t,span=alfa.optim,degree=1,control=loess.control(surface="direct")),data.frame(t=tnuevo),se=FALSE)
Ttnuevo=ts(Ttnuevo,freq=4,start=c(2016,4))#convirtiendo en serie de tiempo al pron?ostico de Tt
ytpron4=Ttnuevo*Stnuevo
ytpron4
accuracy(ytpron4,ytnuevo)
dim(ytpron4)

##### Estabilidad del modelo

##modelo exponencial polinomial estacional

modelo.aux=lm(log(yt)~t+I(t^2)+I(t^3)+I1+I2+I3)
summary(modelo.aux)

##Creando funci?on para estimaci?on recursiva del modelo
estim=function(n){
modelo2=lm(log(Datos)[1:n]~t[1:n]+I(t^2)[1:n]+I(t^3)[1:n]+I1[1:n]+I2[1:n]+I3[1:n]) 
resul1=cbind(coef(modelo2),confint(modelo2))
resul1 
}

#Definicion de tama?o de muestra para las regresiones recursivas
n=matrix(8:length(Datos),ncol=1) 
 

#Realizando las regresiones recursivas. En el arreglo se guardan matrices #cada una de 14 filas (se estiman p=7 par?ametros) y 3 columnas (estimaci?on, LIC, LSC)
b=array(apply(n,1,estim),dim=c(7,3,nrow(n))) 
b

betas0=b[1,,] #Extrayendo informaci?on de los interceptos
betas1=b[2,,] #Extrayendo informaci?on de los betas1
betas2=b[3,,] #Extrayendo informaci?on de los betas2
betas3=b[4,,] #Extrayendo informaci?on de los betas3
Ind_1=b[5,,] #Extrayendo informaci?on de los deltas1
Ind_2=b[6,,] #Extrayendo informaci?on de los deltas2
Ind_3=b[7,,] #Extrayendo informaci?on de los deltas3


#Graficando los beta0 y sus I.C del 95%
win.graph()
matplot(n,t(betas0),type="l",lty=c(1,2,2),col=1,lwd=1.5,ylab=expression(hat(beta)[0])) 
abline(h=coef(modelo.aux)[1],lty=1,col=2,lwd=2)
legend("topright",legend=expression(hat(beta)[0]==1.080e+01),lwd=2,col=2)

#Graficando los beta1 y sus I.C del 95%. Se ajust? rango eje vertical entre 0.005 y 0.011 para ver
win.graph()
matplot(n,t(betas1),type="l",lty=c(1,2,2),ylim=c(-0.1,0.1),col=1,lwd=1.5,ylab=expression(hat(beta)[1]))
abline(h=coef(modelo.aux)[2],lty=1,col=2,lwd=2)
legend("topright",legend=expression(hat(beta)[1]==2.159e-02),lwd=2,col=2)

#Graficando los beta2 y sus I.C del 95%. Se ajust? rango eje vertical entre -3e-05 y 1.e-05 para ver
win.graph()
#con m?s claridad c?mo var?a estimaci?n 
matplot(n,t(betas2),type="l",lty=c(1,2,2),ylim=c(-0.01,0.01),col=1,lwd=1.5,ylab=expression(hat(beta)[2]))
abline(h=coef(modelo.aux)[3],lty=1,col=2,lwd=2)
legend("topright",legend=expression(hat(beta)[2]==1.905e-04),lwd=2,col=2)

#Graficando los beta3 y sus I.C del 95%. Se ajust? rango eje vertical entre -3e-05 y 1.e-05 para ver
win.graph()
matplot(n,t(betas3),type="l",lty=c(1,2,2),col=1,lwd=1.5,ylim=c(-0.00005,0.00005),ylab=expression(hat(beta)[3]))
abline(h=coef(modelo.aux)[4],lty=1,col=2,lwd=2)
legend("topright",legend=expression(hat(beta)[3]==-2.809e-06),lwd=2,col=2)

#Graficando los Ind_1 y sus I.C del 95%
win.graph()
matplot(n,t(Ind_1),type="l",lty=c(1,2,2),col=1,lwd=1.5,ylab=expression(hat(delta)[1]))
abline(h=coef(modelo.aux)[5],lty=1,col=2,lwd=2)
legend("topright",legend=expression(hat(delta)[1]==-1.483e-01),lwd=2,col=2)

#Graficando los Ind_2 y sus I.C del 95%
win.graph()
matplot(n,t(Ind_2),type="l",lty=c(1,2,2),col=1,lwd=1.5,ylab=expression(hat(delta)[2]))
abline(h=coef(modelo.aux)[6],lty=1,col=2,lwd=2)
legend("topright",legend=expression(hat(delta)[2]==-1.023e-01),lwd=2,col=2)

#Graficando los Ind_3 y sus I.C del 95%
win.graph()
matplot(n,t(Ind_3),type="l",lty=c(1,2,2),col=1,lwd=1.5,ylab=expression(hat(delta)[3]))
abline(h=coef(modelo.aux)[7],lty=1,col=2,lwd=2)
legend("topright",legend=expression(hat(delta)[3]==-8.075e-02),lwd=2,col=2)



rr2 =recresid(modelo.aux)
plot(rr2, type = "l",ylab="Residuales recursivos",xlab="t") 
abline(h=0,col=2,lwd=2)
plot(efp(log(yt)~t+I(t^2)+I(t^3)+I1+I2+I3))
sctest(log(yt)~t+I(t^2)+I(t^3)+I1+I2+I3,type="Rec-CUSUM")



plot(ytnuevo,main="Valores reales y pron?sticos ?ltimos Q=4 trimestres\nModelos 1,2,3 y 4",type="b",pch=1,lty=1,col=1,lwd=2,  
    ylab="v?rtice en miles de millones de pesos",ylim=c(min(ytnuevo,pronmod1,pronmod2,ytpron4,ytpronholt),max(ytnuevo,pronmod1,pronmod2,ytpron4,ytpronholt)),xaxt="n")
lines(pronmod1,col=2,pch=2,lty=2,type="b",lwd=2)
lines(pronmod2,col=3,pch=3,lty=3,type="b",lwd=2) 
lines(ytpron4,col=4,pch=4,lty=4,type="b",lwd=2) 
lines(ytpronholt,col=5,pch=5,lty=5,type="b",lwd=2) 
legend("topright",legend=c("Observado","Pron. modelo 1","Pron. modelo 2","Pron. Descomp. & LOESS lineal","Pron. Holt Winters"),
bty="n",col=c(1,2,3,4,5),pch=1:5,lty=c(1,2,3,4,5),lwd=2) 
axis(1,at=time(ytnuevo),labels=c("Trimestre1","Trimestre2","Trimestre3","Trimestre4")) 
 












