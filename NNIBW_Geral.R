library(ggplot2)
library(survival)
library(reshape2)
library(dplyr)
library(gridExtra)
theme_set(theme_bw())
getwd()

# FUNÇÕES DO MODELO -----

# fd  sem o m:
fd_WB <- function(x,beta=2,lambda=1,delta=0.5,mu=-1){
  
  if(mu>=0){
    cond <- mu^((1)/(delta+1))
    f <- ifelse(x>=cond,(beta*(delta+1)/lambda)*exp(-((x^(delta+1)-mu)/lambda)^beta)*((x^(delta+1)-mu)/lambda)^(beta-1)*x^delta,
                0)    
  }
  else{
    cond <- -(-mu)^(1/(delta+1))
    a1 <- ((beta*(delta+1))/lambda)
    a2 <- exp(-((x^(delta+1)-mu)/lambda)^beta)
    a3 <- (((x)^(delta+1)-mu)/lambda)^(beta-1)
    a = a1*a2*a3*x^delta
    
    b1 <- (beta*(delta+1)/lambda)
    b2 <- exp(-((-(-x)^(delta+1)-mu)/lambda)^beta)
    b3 <-  ((-(-x)^(delta+1)-mu)/lambda)^(beta-1)
    b = b1*b2*b3*(-x)^delta
    
    f <- ifelse(x>=0,a,
                ifelse(x<0 & x>cond,b,
                       0))
  }
  
  # print(c(a1,a2,a3,b1,b2,b3))
  return(f)
}


# Função adicionando um novo parâmetro m: (1.2.4)
fd_WBm <- function(y,beta,lambda,delta,mu,m){
  #m =  (-mu)^(1/(delta+1))
  f <- fd_WB(x=(y-m),beta,lambda,delta,mu)
  return(f)
}

# Função caso particular (1.2.22)
fd_WBm_part <- function(y,beta,lambda,delta,mu){
  m =  (-mu)^(1/(delta+1))
  f <- fd_WB(x=(y-m),beta,lambda,delta,mu)
  return(f)
}


# fdA original sem o m (Acumulada):
FD_WB <- function(x,beta=2,lambda=1,delta=0.5,mu=-1){
  
  if(mu>=0){
    cond <- mu^((1)/(delta+1))
    a2 <- exp(-((x^(delta+1)-mu)/lambda)^beta)
    a = 1-a2
    f <- ifelse(x>=cond,a,
                0)    
  }
  else{
    cond <- -(-mu)^(1/(delta+1))
    a2 <- exp(-((x^(delta+1)-mu)/lambda)^beta)
    a = 1-a2
    
    b2 <- exp(-((-(-x)^(delta+1)-mu)/lambda)^beta)
    b = 1-b2
    
    f <- ifelse(x>=0,a,
                ifelse(x<0 & x>cond,b,
                       0))
  }
  
  return(f)
}


# Função adicionando um novo parâmetro m + (caso particular ou não), (1.2.5) e (1.2.21):
FD_WBm <- function(y,beta,lambda,delta,mu){
  m =  (-mu)^(1/(delta+1)) # comentar essa linha se não for o caso particular
  f <- FD_WB(x=(y-m),beta,lambda,delta,mu)
  
  return(f)
}

# função quantil  caso particular:

fquantil <- function(y,beta=1,lambda=1,delta=0,mu=0){
  m = (-mu)^(1/(1+delta)) 
  cond <- (mu + lambda*((-log(1-y))^(1/beta)))
  a <- (mu + lambda*((-log(1-y))^(1/beta)))^(1/(delta+1)) -m
  b <- -(-mu -lambda*((-log(1-y))^(1/beta)))^(1/(delta+1)) -m 
  
  f <- ifelse(cond>=0,a,b)
  return(f)
}

fquantilm_part <- function(y,beta=1,lambda=1,delta=0,mu=0){
  cond <- (mu + lambda*((-log(1-y))^(1/beta)))
  m=  (-mu)^(1/(delta+1))
  a <- (mu + lambda*((-log(1-y))^(1/beta)))^(1/(delta+1))+m
  b <- -(-mu-lambda*((-log(1-y))^(1/beta)))^(1/(delta+1))+m
  
  f <- ifelse(cond>=0,a,ifelse(cond<0 ,b,0))
  
  return(f)
}

# funcao tempo de retorno

zpfun <- function(p,beta=1,lambda=1,delta=0,mu=0){
  m =  (-mu)^(1/(delta+1))
  cond <- exp(-(-(mu/lambda))^beta)
  a <- (mu + lambda*((-log(p))^(1/beta)))^(1/(delta+1)) + m
  b <- -(-mu -lambda*((-log(p))^(1/beta)))^(1/(delta+1)) + m
  
  f <- ifelse(p<=cond,a,b)
  return(f)
}



# função de risco caso particular 

h_WB_part <- function(y,beta=2,lambda=1,delta=0.5,mu=-1){
  m =  (-mu)^(1/(delta+1))
  x = (y-m)
  
  
  cond <- -(-mu)^(1/(delta+1))
  a1 <- ((beta*(delta+1))/lambda)
  a2 <- 1
  a3 <- (((x)^(delta+1)-mu)/lambda)^(beta-1)
  a = a1*a2*a3*x^delta
  
  b1 <- ((beta*(delta+1))/lambda)
  b2 <- 1
  b3 <-  ((-(-x)^(delta+1)-mu)/lambda)^(beta-1)
  b = b1*b2*b3*(-x)^delta
  
  f <- ifelse(x>=0,a,
              ifelse(x<0 & x>=cond,b,
                     0))
  
  return(f)
}


# Log verossimilhança:

# log-vero geral
lvero_geral <- function(par){
  beta <-par[1]
  lambda <-par[2]
  delta<-par[3]
  mu <-par[4]
  m <- par[5]
  z = x - m
  n <- length(x)
  cond1 <- mu^(1/(1+delta))
  cond2 <- -(-mu)^(1/(1+delta))
  
  if(mu>=0){
    # mu>0
    #L1 (sem a soma e o n, pois soma no final)
    b1 <- log(beta/lambda)
    b2 <- - ((z^(delta+1)-mu)/lambda)^beta
    b3 <- (beta-1)*(log((z^(delta+1)-mu)/lambda))
    b4 <- log((delta+1)*z^delta)
    b = b1+b2+b3+b4 
    
    f <- ifelse(z<=cond1,0,b)
    
  }
  
  
  else{ #mu<0
    # L2
    a1 <- log(beta/lambda)
    a2 <- - ((-(-z)^(delta+1)-mu)/lambda)^beta
    a3 <- (beta-1)*(log((-(-z)^(delta+1)-mu)/lambda))
    a4 <- log((delta+1)*(-z)^delta)
    a = a1+a2+a3+a4
    
    b1 <- log(beta/lambda)
    b2 <- - ((z^(delta+1)-mu)/lambda)^beta
    b3 <- (beta-1)*(log((z^(delta+1)-mu)/lambda))
    b4 <- log((delta+1)*z^delta)
    b = b1+b2+b3+b4
    
    # print('as');print(a2);print(a3);print(a4);
    
    # print('bs');print(b2);print(b3);print(b4);
    
    f <- ifelse(z<=cond2,0,
                ifelse(z>cond2 & z<0,a,
                       b))
    
  }
  #print(f)
  #print(log(w1));print(log(w2))
  return(f)#sum(f))
}


# caso particular x estrit positivo m= mu...

lvero_part <- function(par){
  beta <-par[1]
  lambda <-par[2]
  delta<-par[3]
  mu <-par[4]
  m <- (-mu)^(1/(delta+1))
  v = x - m
  
  # L1
  ind1 <- ifelse(x<m,1,0)

  #L2
  ind2 <- ifelse(x>m,1,0)
  
  
  # L1
  options(warn=-1) 
  a1 <- log(beta/lambda) +log(delta+1)
  a2 <- - ((-(-v)^(delta+1)-mu)/lambda)^beta
  a3 <- (beta-1)*(log((-(-v)^(delta+1)-mu)/lambda))
  a4 <- delta*log(-v)
  a = a1+a2+a3+a4
  
  #L2
  b1 <- log(beta/lambda) +log(delta+1) 
  b2 <- - ((v^(delta+1)-mu)/lambda)^beta
  b3 <- (beta-1)*(log((v^(delta+1)-mu)/lambda))
  b4 <- delta*log(v)
  b = b1+b2+b3+b4
  
  # Coloca 0 no NaN (talvez nem precisasse disso,
  # só remover e somar)
  b[is.nan(b)] <- 0   
  a[is.nan(a)] <- 0  
  b[(b==-Inf)] <- 0   
  a[a==-Inf] <- 0  
  
  
  f = sum(ind1*a,ind2*b)
  
  
  # print(a)
  # print(b)
  # print(ind1)
  # print(ind2)
  return(f)#sum(f))
}

# condicionando o limite dos parametros
lvero_part_2 <- function(par){
  beta <-par[1]
  lambda <-par[2]
  delta<-par[3]
  mu <-par[4]
  m <- (-mu)^(1/(delta+1))
  v = x - m
  
  if((beta>0) && (lambda>0) && (delta>0) && (mu<0) ) return((-1)*lvero_part(c(beta,lambda,delta,mu)))
  
  else return(-Inf)
  
}

# ----------------- Função do prof. Vila
weibivar <- function(x,theta){
  alpha <- theta[1]
  beta <- theta[2]
  delta <- theta[3]
  
  zt = 2 + (delta^2)*(beta^2)*gamma(1+(2/alpha)) - 2*delta*beta*gamma(1+(1/alpha))
  
  a1 = (alpha/(beta*zt))*(1 + (1 - delta*x)^2)
  a2 = (x/beta)^(alpha-1)*exp(-(x/beta)^alpha)
  
  return(a1*a2)
}

##### Função de Verossimilhança do Vila
vero<-function(theta,x){
  alpha<-theta[1]
  beta<-theta[2]
  delta<-theta[3]
  if ((alpha>0)&&(beta>0)) return (-1*sum(log(weibivar(x,theta))))
  else return (-Inf)
}



# EXEMPLO PLOT FUNÇÃO ----

t <- seq(-3,3,0.05)

# alterando o m  # ORDEM beta lambda delta mu m
valor=fd_WBm(t,2,2,2,-1,-1)
valor.1=fd_WBm(t,2,2,2,-1,-2)
valor.2=fd_WBm(t,2,2,2,-1,0)
valor.3=fd_WBm(t,2,2,2,-1, 1)

# beta lmabda delta mu
plotar <- data.frame(valor,t)
plotar2 <- data.frame(valor.1,t)
plotar3 <- data.frame(valor.2,t) ##
plotar4 <- data.frame(valor.3,t)
plotar <- data.frame(plotar,plotar2,plotar3,plotar4)

library(ggplot2)
theme_set(theme_bw())

p1 <- ggplot(plotar,aes(t,valor)) + 
  geom_line(aes(color=" m=-1"),size=1) +
  
  geom_line(aes(t.1,valor.1,color=" m=-2"),size=1) +
  
  geom_line(aes(t.2,valor.2,color=" m=0"),size=1) +
  
  geom_line(aes(t.3,valor.3,color=" m=1"),size=1) +
  
  labs(color="parâmetros",
       y="f(x)",x='x') + theme(legend.position = c(.80, .75),
                               legend.box.background = element_rect(color="black", size=1),
                               legend.key.size = unit(0.2, 'cm'))+
  xlim(-3,6) + 
  theme(legend.title = element_blank())

p1


# SIMULAÇÃO ----


b=c(1,2);l=c(1,2);d=c(.5,1);m=c(-1,-2)

df <- expand.grid(b,l,d,m)

df <- df[-c(10,14),]

n <- 500 #tamanho da amostra
rep <- 1000 #replicacoes montecarlo

resultpar <- matrix(,nrow(df),4)
eqm <- matrix(,nrow(df),4)

for(i in 1:nrow(df)){
  beta <- df[i,1]
  lambda <- df[i,2]
  delta <- df[i,3]
  mu <- df[i,4]
  
  mat_temp <- matrix(,rep,4)
  eqm_temp <- matrix(,rep,4)
  
  for(j in 1:rep){
    
    y <- runif(n)
    x <- fquantilm_part(y,beta=beta,lambda=lambda,delta=delta,mu=mu) 
    
    tryCatch(
      error = function(cnd) NA,
      ML_est2 <-  optim(c(beta,lambda,delta,mu) , lvero_part,
                        control=list(fnscale=-1,maxit=1000),
                        lower=c(0.0001,0.0001,0.0001,-Inf), 
                        upper=c(rep(Inf,3),0.0001),
                        method="L-BFGS-B",
                        hessian=T)
    )
    
    
    
    mat_temp[j,] <- ML_est2$par
    eqm_temp[j,] <- as.numeric((ML_est2$par - df[i,])^2)
  }
  
  resultpar[i,] <- colMeans(mat_temp)
  eqm[i,] <- colMeans(eqm_temp)
  
}

c(t(round(resultpar,3)))
c(t(round(eqm,3)))
#format(round(eqm,3), scientific = FALSE)

# Comentar e descomentar aqui dependendo do tamanho da amostra que estão sendo geradas para o gráfico
# resultpar50 <- resultpar
# eqm50 <- eqm

resultpar500 <- resultpar
eqm500 <- eqm

library(ggplot2)
theme_set(theme_bw())

# Gráficos do apêndice 
for(i in 1:nrow(df)){
  x <- seq(0,5,0.05)
  original = fd_WBm_part(x,df[i,1],df[i,2],df[i,3],df[i,4])
  estimado50 = fd_WBm_part(x,resultpar50[i,1],resultpar50[i,2],resultpar50[i,3],resultpar50[i,4])
  estimado500 = fd_WBm_part(x,resultpar500[i,1],resultpar500[i,2],resultpar500[i,3],resultpar500[i,4])
  
  plotar <- data.frame(original,estimado50,estimado500,x)
  
  p1 <- ggplot(plotar,aes(x,original)) + 
    geom_line(aes(color='Theorical')) +
    geom_point(aes(x,estimado50,color='Fitted: n = 50'),size=.7) +
    geom_point(aes(x,estimado500,color='Fitted: n = 500'),size=.7) +
    labs(y="f(x)",x='x') + theme(legend.position = c(.70, .75),
                                 legend.box.background = element_rect(color="black", size=1),
                                 legend.key.size = unit(0.1, 'cm')) +
    scale_color_manual(values=c("cornflowerblue", "firebrick1","gray7")) + 
    theme(legend.title = element_blank())
  #p1
  
  ggsave(filename = paste("imagensR\\simulacaocurvascompar","_",i,".png",sep=""),p1,
         width = 5, height = 3, dpi = 300, units = "in", device='png')
  
}



# ESTIMAÇÃO DADOS DELHI ----

novadelhi <- read.csv("dados\\DailyDelhiClimateTest.csv",encoding="UTF-8",header=T)
#
head(novadelhi)
hist(novadelhi$meantemp, breaks = 12)

x <- novadelhi$meantemp
x <- exp(x/100) 
hist(x)
x <- x[!is.na(x)]
hist(x)
anyNA(x)

df <- cbind(runif(1000,8,12),runif(1000,.5,2),runif(1000,.05,1),-runif(1000,0.5,1.5))

valormax <- c()
estimados <- list()
vari <- list()
codigo <- c()
confirm <- c()

# estimação pra cada chute inicial
for(i in 1:nrow(df)){
  tryCatch({
    a <- optim(df[i,] , lvero_part,
               control=list(fnscale=-1,maxit=10000),
               lower=c(0.0001,0.0001,0.03,-Inf), 
               upper=c(rep(Inf,3),-0.01),
               method="L-BFGS-B",
               hessian=T) 
    valormax[i] <- a$value
    estimados[[i]] <- a$par
    vari[[i]] <- solve((-1)*a$hessian)
    confirm[i] <- any(diag(solve((-1)*a$hessian))<0)
    codigo[i] <- a$convergence
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  cat(i,a$value,"\n")
}

# valores de máximo estimados: 
head(sort(valormax,decreasing=T),10)

dest <- head(estimados[order(valormax,decreasing = T)],10)
quemsao <- head(df[order(valormax,decreasing = T),],10)

# valores dos parâmetros nos máximos:
dest

# verificando a convergencia e se a matriz hessiana é negativa:
head(codigo[order(valormax,decreasing = T)],10)
head(confirm[order(valormax,decreasing = T)],10)

# ESTIMAÇÃO DADOS SVALBARD ----


teste <- read.csv("dados/svalbard-climate-1912-2017.csv",encoding="UTF-8",header=T)#,skip=9)
#
head(teste)

library(reshape2)

#organizando os dados:
teste <- melt(teste,id = "YEAR")
svalbar <- teste[teste$value!=999.9,]

hist(svalbar$value)
x <- (svalbar$value)


# transformação pra deixar positiva
x <- exp(x/10)

anyNA(x)
hist(x)


df <- cbind(runif(1000,8,12),runif(1000,.5,2),runif(1000,.05,1),-runif(1000,0.5,1.5))

valormax <- c()
estimados <- list()
vari <- list()
codigo <- c()
confirm <- c()

for(i in 1:nrow(df)){
  tryCatch({
    a <- optim(df[i,] , lvero_part,
               control=list(fnscale=-1,maxit=10000),
               lower=c(0.0001,0.0001,0.03,-Inf), 
               upper=c(rep(Inf,3),-0.01),
               method="L-BFGS-B",
               hessian=T) 
    valormax[i] <- a$value
    estimados[[i]] <- a$par
    vari[[i]] <- solve((-1)*a$hessian)
    confirm[i] <- any(diag(solve((-1)*a$hessian))<0)
    codigo[i] <- a$convergence
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  cat(i,a$value,"\n")
}


# valores de máximo estimados: 
head(sort(valormax,decreasing=T),10)

dest <- head(estimados[order(valormax,decreasing = T)],10)
quemsao <- head(df[order(valormax,decreasing = T),],10)

# valores dos parâmetros nos máximos:
dest

# verificando a convergencia e se a matriz hessiana é negativa:
head(codigo[order(valormax,decreasing = T)],10)
head(confirm[order(valormax,decreasing = T)],10)


# ESTIMAÇÃO DADOS YELLOWKNIFE  ----

# https://climate.weather.gc.ca/

arquivos <- list.files(path = "dados\\yellowknifeanos",full.names=T)

yknife <- read.csv(arquivos[1],encoding="UTF-8",header=T)
for(i in 2:length(arquivos)){
  a <- read.csv(arquivos[i],encoding="UTF-8",header=T)
  yknife <- rbind(yknife,a)
  
}

head(yknife)
tail(yknife)

yknife$Season <- factor(yknife$Month, levels = c(1,2,3,4,5,6,7,8,9,10,11,12), 
                        labels = c("Winter","Winter","Winter", "Spring", "Spring", "Spring",
                                   "Summer","Summer","Summer","Fall","Fall","Fall"))



x <- yknife$Min.Temp...C.
hist(x)
x <- exp(x/100) 
yknife <- yknife[!is.na(x),]
x <- x[!is.na(x)]
seas <- yknife$Season

N<-length(x)
n<-90
tau <- floor(N/n)
M<-matrix(0,tau,1)
Mseason <- c()
j<-1
for (i in 1:tau){
  M[i]<-min(x[j:(j+n-1)])
  Mseason[i] <- seas[j:(j+n-1)][which.min(x[j:(j+n-1)])]
  j<-j+n}
head(M)
sum(is.na(M))
hist(M, prob=T)
lines(density(M),lwd=2)
x <- M


hist(x)

df <- cbind(runif(1000,2,30),runif(1000,.05,40),runif(1000,.001,10),-runif(1000,0.02,30))

valormax <- c()
estimados <- list()
vari <- list()
codigo <- c()
confirm <- c()

for(i in 1:nrow(df)){
  tryCatch({
    a <- optim(df[i,] , lvero_part,
               control=list(fnscale=-1,maxit=10000),
               lower=c(0.0001,0.0001,0.03,-Inf), 
               upper=c(rep(Inf,3),-0.01),
               method="L-BFGS-B",
               hessian=T) 
    valormax[i] <- a$value
    estimados[[i]] <- a$par
    vari[[i]] <- solve((-1)*a$hessian)
    confirm[i] <- any(diag(solve((-1)*a$hessian))<0)
    codigo[i] <- a$convergence
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  cat(i,a$value,"\n")
}


# valores de máximo estimados: 
head(sort(valormax,decreasing=T),10)

dest <- head(estimados[order(valormax,decreasing = T)],10)
quemsao <- head(df[order(valormax,decreasing = T),],10)

# valores dos parâmetros nos máximos:
dest

# verificando a convergencia e se a matriz hessiana é negativa:
head(codigo[order(valormax,decreasing = T)],10)
head(confirm[order(valormax,decreasing = T)],10)



# ESTIMAÇÃO DADOS DAWSON CITY ----


dcity <- read.csv("dados/en_climate_daily_YT_2100407_2020_P1D.csv",encoding="UTF-8",header=T)#,skip=9)
#
head(dcity)
names(dcity)

hist(dcity$Max.Temp...C., breaks = 12)
hist(dcity$Spd.of.Max.Gust..km.h., breaks = 12)

x <- dcity$Max.Temp...C.
x <- exp(x/100) # muda ser sobre 10 ou 100
x <- x[!is.na(x)]
# x <- sample(x,100)
hist(x)

df <- cbind(runif(1000,1,5),runif(1000,.05,4),
            runif(1000,.05,4),-runif(1000,0.1,3))

valormax <- c()
estimados <- list()
vari <- list()
codigo <- c()
confirm <- c()


for(i in 1:nrow(df)){
  tryCatch({
    a <- optim(df[i,] , lvero_part,
               control=list(fnscale=-1,maxit=10000),
               lower=c(0.0001,0.0001,0.03,-Inf), 
               upper=c(rep(Inf,3),-0.01),
               method="L-BFGS-B",
               hessian=T) 
    valormax[i] <- a$value
    estimados[[i]] <- a$par
    vari[[i]] <- solve((-1)*a$hessian)
    confirm[i] <- any(diag(solve((-1)*a$hessian))<0)
    codigo[i] <- a$convergence
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  cat(i,"\n")
}

# valores de máximo estimados: 
head(sort(valormax,decreasing=T),10)

dest <- head(estimados[order(valormax,decreasing = T)],10)
quemsao <- head(df[order(valormax,decreasing = T),],10)

# valores dos parâmetros nos máximos:
dest

# verificando a convergencia e se a matriz hessiana é negativa:
head(codigo[order(valormax,decreasing = T)],10)
head(confirm[order(valormax,decreasing = T)],10)



# PLOT DADOS GERAIS ----

# --

s1 <- ggplot(svalbar, aes(value)) +
  geom_histogram(aes(y = ..density..),fill='grey',col='grey4',binwidth = 2) +
  labs(subtitle='Svalbard',x='Temperature C°')


s2 <- ggplot(yknife2, aes(maxx)) +
  geom_histogram(aes(y = ..density..),fill='grey',col='grey4',bins = 10) +
  labs(subtitle='Yellow Knife',x='Temperature C°')


s3 <- ggplot(dcity, aes(Max.Temp...C.)) +
  geom_histogram(aes(y = ..density..),fill='grey',col='grey4',binwidth = 5) +
  labs(subtitle='Dawson City',x='Temperature C°')

s4 <- ggplot(novadelhi, aes(meantemp)) +
  geom_histogram(aes(y = ..density..),fill='grey',col='grey4',binwidth = 1.5) +
  labs(subtitle='Delhi',x='Temperature C°')

require(gridExtra)
grid.arrange(s1, s2,s3,s4, ncol=2)

ggsave("imagensR/hist_originais0.png", arrangeGrob(s1, s2,s3,s4),
       width = 5, height = 3, dpi = 300, units = "in", device='png')




# PLOT DADOS AJUSTADOS ----

#dcity
x <- exp(dcity$Max.Temp...C./100)
x <- x[!is.na(x)]
dcity <- dcity[!is.na(dcity$Max.Temp...C.),]

beta=  19.594493726281
lambda=  1.10251321448119
delta =  0.675155986643908
mu=  -1.08832025976928
dcity$y = fd_WBm_part(x, beta = beta,lambda=lambda, delta =delta, mu = mu)


# estimação Vila
chute<-c(2.5,2.5,0.2)

nlmmodel <- nlm(vero,chute,x=x)
nlmmodel

theta_est <- nlmmodel$estimate
dcity$prof=weibivar(x,theta_est)



ggplot(dcity, aes(exp(Max.Temp...C./100))) +
  geom_histogram(aes(y = ..density..),fill='grey',col='grey4',binwidth = .05) +
  geom_line(aes(exp(Max.Temp...C./100),y),col="blue2") +
  geom_line(aes(exp(Max.Temp...C./100),prof),col="red") +
  labs(subtitle='Dawson City',x='y')



#yknife 
x <- exp(yknife$Min.Temp...C./100)
x <- x[!is.na(x)]
yknife <- yknife[!is.na(yknife$Min.Temp...C.),]


N<-length(x)
n<-90
tau <- floor(N/n)
M<-matrix(0,tau,1)
j<-1
for (i in 1:tau){
  M[i]<-min(x[j:(j+n-1)])
  j<-j+n}
head(M)
sum(is.na(M))
hist(M, prob=T)
lines(density(M),lwd=2)
x <- M

summary(log(M)*100)


#minimos

beta = 84.2695289
lambda =  0.6606731
delta =1.6552931
mu = -0.6626890


yknife2 = data.frame(x)
yknife2$maxx = 100*log(x)
yknife2$y = fd_WBm_part(x, beta = beta,lambda=lambda, delta =delta, mu = mu)

chute<-c(2.5,2.5,0.2)

nlmmodel <- nlm(vero,chute,x=x)
nlmmodel

theta_est <- nlmmodel$estimate
yknife2$prof=weibivar(x,theta_est)


z2 <- ggplot(yknife2, aes(x)) +
  geom_histogram(aes(y = ..density..),fill='grey',col='grey4',bins = 9) +
  geom_line(aes(x,y),col="blue2") +
  geom_line(aes(x,prof),col="red") +
  labs(subtitle='Yellow Knife',x='y')

z2


#sva
x <- exp(svalbar$value/10)
x <- x[!is.na(x)]
beta=  2.09820224047842
lambda=  1.02254548476456
delta =   0.541382456353869
mu=  -1.10301957817666

svalbar$y = fd_WBm_part(x, beta = beta,lambda=lambda, delta =delta, mu = mu)
chute<-c(1.5,2,0.7)

nlmmodel <- nlm(vero,chute,x=x)
theta_est <- nlmmodel$estimate
svalbar$prof=weibivar(x,c(1.5,1.5,1))

ggplot(svalbar, aes(exp(value/10))) +
  geom_histogram(aes(y = ..density..),fill='grey',col='grey4',binwidth = .15) +
  geom_line(aes(exp(value/10),y),col="blue2") +
  geom_line(aes(exp(value/10),prof),col="red") +
  labs(subtitle='Svalbard',x='y')



#delhi
x <- exp(novadelhi$meantemp/100)

beta=  201.897361761866
lambda=  1.68789121065167
delta =   1.01582450161826
mu=  -1.68984633468315

novadelhi$y = fd_WBm_part(x, beta = beta,lambda=lambda, delta =delta, mu = mu)
chute<-c(1,1.5,2)

nlmmodel <- nlm(vero,chute,x=x)
theta_est <- nlmmodel$estimate
novadelhi$prof=weibivar(x,theta_est)


# todos os plot:

z1 <- ggplot(svalbar, aes(exp(value/10))) +
  geom_histogram(aes(y = ..density..),fill='grey',col='grey4',binwidth = .15) +
  geom_line(aes(exp(value/10),y),col="blue2") +
  geom_line(aes(exp(value/10),prof),col="red") +
  labs(subtitle='Svalbard',x='y')


z2 <- ggplot(yknife2, aes(exp(maxx/100))) +
  geom_histogram(aes(y = ..density..),fill='grey',col='grey4',bins = 9) +
  geom_line(aes(exp(maxx/100),y),col="blue2") +
  geom_line(aes(exp(maxx/100),prof),col="red") +
  labs(subtitle='Yellow Knife',x='y')



z3 <- ggplot(dcity, aes(exp(Max.Temp...C./100))) +
  geom_histogram(aes(y = ..density..),fill='grey',col='grey4',binwidth = .05) +
  geom_line(aes(exp(Max.Temp...C./100),y),col="blue2") +
  geom_line(aes(exp(Max.Temp...C./100),prof),col="red") +
  labs(subtitle='Dawson City',x='y')


z4 <- ggplot(novadelhi, aes(exp(meantemp/100))) +
  geom_histogram(aes(y = ..density..),fill='grey',col='grey4',binwidth = .03) +
  geom_line(aes(exp(meantemp/100),y),col="blue2") +
  geom_line(aes(exp(meantemp/100),prof),col="red") +
  labs(subtitle='Delhi',x='y')




ggsave("imagensR/hist_transform.png", arrangeGrob(z1,z2,z3,z4),
       width = 5, height = 3, dpi = 300, units = "in", device='png')


# PLOT DADOS INDIVIDUAIS ----


# dawson city -----------------------
getwd()
# C:\\Users\\Cris\\Documents\\Bia\\Mono\\R
# D:\\foque\\unb\\!Mestrado\\Mono\\R\\

dcity <- read.csv("dados\\en_climate_daily_YT_2100407_2020_P1D.csv",encoding="UTF-8",header=T)#,skip=9)
#
head(dcity)
x <- dcity$Max.Temp...C.
x <- exp(x/100) 
x <- x[!is.na(x)]
hist(x)
max( dcity$Max.Temp...C.,na.rm=T)
sort(( dcity$Max.Temp...C.))
beta=  19.594493726281
lambda=  1.10251321448119
delta =  0.675155986643908
mu=  -1.08832025976928


# LOGVERO 163.637
#AIC: AIC = - 2*log L + k * GL, k=2

-2*163.637 + 2*4

# modelo prof Vila:
chute<-c(1.5,2,0.7)
nlmmodel <- nlm(vero,chute,x=x)
theta_est <- nlmmodel$estimate
prof=weibivar(x,theta_est)

-2*((-1)*nlmmodel$minimum) + 2*4



df <- data.frame(x = x,y = fd_WBm_part(x, beta = beta,lambda=lambda, delta =delta, mu = mu),
                 prof=weibivar(x,theta_est),
                 date=c(1:length(x)))#dcity$Date.Time[!is.na(dcity$Max.Temp...C.)])

head(df)

g1 <- ggplot(df, aes(x)) + geom_histogram(aes(y = ..density..),fill="gray75", 
                                          binwidth = .05, 
                                          col="black", 
                                          size=.1) +
  geom_line(aes(x,y),col="blue4") +labs(x="y") 

g1

# serie temporal de x
g2 <- ggplot(df) +
  geom_line(aes(date,x),col="blue4") + scale_x_continuous(breaks=cumsum(c(31,28,31,30,31,30,31,31,30,31,30,31)),
                                                          labels=c("jan","","","","may","","",
                                                                   "","sep","","","dec")) +labs(y="y")
g2



# Acumuladas
sobj <- survfit(Surv(x,rep(1,length(x)))~1)

dsurv <- data.frame(sobj$time,empirical = 1 - sobj$surv,sobj.surv = sobj$surv)

df2 <- merge(df,dsurv,by.x = "x",by.y = "sobj.time")

df2$index <- sort(df2$x)

df2$theorical <- FD_WBm(sort(df$x), beta = beta,lambda=lambda, delta =delta, mu = mu)

plotar <- df2 %>% select(index,theorical,empirical) %>% melt(id = "index")
tail(plotar)

g3 <- ggplot(df2) +
  geom_line(aes(index,theorical,col='theorical')) + 
  geom_path(aes(index,empirical,col='empirical')) +
  theme(legend.position = c(.20, .85),
        legend.key = element_blank(),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_blank()) +
  scale_colour_manual(values=c(theorical="blue4",empirical="firebrick3")) +
  labs(x="y",y="ECDF")  + guides(color=guide_legend(override.aes=list(fill=NA)))

g3

head(df2)

xxx <- df2$x[-nrow(df2)]
yyy <- fquantilm_part(1-df2$sobj.surv,beta=beta,lambda=lambda,delta=delta,mu=mu)[-nrow(df2)]
plot(xxx,yyy)
points(c(0,3),c(0,3),type='l')


df2$quantile <- sort(fquantilm_part(1-df2$sobj.surv,beta = beta,lambda=lambda, delta =delta, mu = mu))

g4 <- ggplot(df2[-nrow(df2),]) +
  geom_point(aes(index,quantile)) + 
  geom_abline(col='red') +labs(x="y")
g4


require(gridExtra)
grid.arrange(g1, g2,g3,g4, ncol=2)


ggsave("imagensR/grid_dawsoncity.png", arrangeGrob(g2,g1,g3,g4),
       width = 5, height = 3, dpi = 300, units = "in", device='png')



head(df2)

# tempo de retorno n ----
xp <- seq(.3,.9999,0.0001)
TR <- 1/(1-xp)
zpp <- sort(fquantilm_part(xp,beta = beta,lambda=lambda, delta =delta, mu = mu))
# teorico
zpporiginal <- 100*log(zpp)

# empirico
zppemp <- sort(fquantilm_part(1-df2$sobj.surv,beta = beta,lambda=lambda, delta =delta, mu = mu))
zppemporiginal <- 100*log(zppemp)
xpemp <- 1-df2$sobj.surv
TRemp <- 1/(1-xpemp)



# x temperatura e y em dias, colocar legenda emp e theo, plot original e sep
plot(zpp,TR,type="l",ylim = c(0,2000),xlab="Return level X ",ylab="Return time (days)",main="Dawson City")
points(zppemp,pch=16,TRemp,col="red")
legend(1.16, 1600, legend=c("Empirical", "Theorical"),pch = c(19, NA), lty = c(NA, 1),
       col=c("red", "black"), cex=0.8)


plot(zpporiginal,TR,type="l",ylim = c(0,2000),xlim=c(20,31),xlab="Return level (C°)",ylab="Return time (days)",main="Dawson City")
points(zppemporiginal,pch=16,TRemp,col="red")
legend(20, 1600, legend=c("Empirical", "Theorical"),pch = c(19, NA), lty = c(NA, 1),
       col=c("red", "black"), cex=0.8)

which(xp==.375) #quantil # 2001
TR[751] # tempo de retorno
zpporiginal[751] # temperatura

which(xp==.50) #quantil # 2001
TR[2001] # tempo de retorno
zpporiginal[2001] # temperatura


which(xp==.9000) #quantil # 2001
TR[2001] # tempo de retorno
zpporiginal[2001] # temperatura


which(xp==.95) #quantil
TR[which(xp==.95)] # tempo de retorno
zpporiginal[2501] # temperatura

which(xp==.99) #quantil
TR[which(xp==.99)] # tempo de retorno
zpporiginal[2901] # temperatura

which(xp==.9991) #quantil
TR[which(xp==.9991)-1] # tempo de retorno
zpporiginal[2991] # temperatura

which(xp==.9994) #quantil
TR[which(xp==.9994)+1] # tempo de retorno
zpporiginal[2995] # temperatura


# YELLOWKNIFE ----

library(ggplot2)


arquivos <- list.files(path = "dados\\yellowknifeanos",full.names=T)

yknife <- read.csv(arquivos[1],encoding="UTF-8",header=T)
for(i in 2:length(arquivos)){
  a <- read.csv(arquivos[i],encoding="UTF-8",header=T)
  yknife <- rbind(yknife,a)
  
}

head(yknife)
tail(yknife)

yknife$Season <- factor(yknife$Month, levels = c(1,2,3,4,5,6,7,8,9,10,11,12), 
                        labels = c("Winter","Winter","Winter", "Spring", "Spring", "Spring",
                                   "Summer","Summer","Summer","Fall","Fall","Fall"))



x <- yknife$Min.Temp...C.
hist(x)
x <- exp(x/100) # muda ser sobre 10 ou 100
yknife <- yknife[!is.na(x),]
x <- x[!is.na(x)]
seas <- yknife$Season

N<-length(x)
n<-90
tau <- floor(N/n)
M<-matrix(0,tau,1)
Mseason <- c()
j<-1
for (i in 1:tau){
  M[i]<-min(x[j:(j+n-1)])
  Mseason[i] <- seas[j:(j+n-1)][which.min(x[j:(j+n-1)])]
  j<-j+n}
head(M)
sum(is.na(M))
hist(M, prob=T)
lines(density(M),lwd=2)
x <- M


beta =  84.2695289  
lambda =  0.6606731  
delta = 1.6552931 
mu = -0.6626890

yknife2 = data.frame(x)
yknife2$maxx = 100*log(x)
yknife2$y = fd_WBm_part(x, beta = beta,lambda=lambda, delta =delta, mu = mu)

chute<-c(1.5,2.5,0.7)

nlmmodel <- nlm(vero,chute,x=x)
nlmmodel

theta_est <- nlmmodel$estimate
yknife2$prof=weibivar(x,theta_est)


acf(M,main="")

plot(M,type="l")

dest
# LOGVERO 259.0536
#AIC: AIC = - 2*log L + k * GL, k=2

-2*259.0536 + 2*4

#prof vila
chute<-c(1.5,2,0.7)
nlmmodel <- nlm(vero,chute,x=x)
theta_est <- nlmmodel$estimate

-2*((-1)*nlmmodel$minimum) + 2*4

# -295.8803

df <- data.frame(x = x,y = fd_WBm_part(x, beta = beta,lambda=lambda, delta =delta, mu = mu),
                 prof=weibivar(x,theta_est),
                 date=c(1:length(x)))

head(df)

g1 <- ggplot(df, aes(x)) + geom_histogram(aes(y = ..density..),fill="gray75", 
                                          binwidth = .05, 
                                          col="black", 
                                          size=.1) +
  geom_line(aes(x,y),col="blue4") +labs(x="y") 

g1

38/4

anos <- rep(c(1953,2014,2015,2016,2017,2018,2019,2020,2021,2022),each=4)
df$anos <- anos[-c(39,40)]
# serie temporal de x
g2 <- ggplot(df) +
  geom_line(aes(date,x),col="blue4") + scale_x_continuous(breaks=(c(2,100,200,279)),
                                                          labels=c("1953","1978","2003","2022")) +labs(y="y")
g2 

nrow(df)

library(survival)
library(dplyr)
library(reshape2)
# Acumuladas
sobj <- survfit(Surv(x,rep(1,length(x)))~1)

dsurv <- data.frame(sobj$time,empirical = 1 - sobj$surv,sobj.surv = sobj$surv)

df2 <- merge(df,dsurv,by.x = "x",by.y = "sobj.time")

df2$index <- sort(df2$x)

df2$theorical <- FD_WBm(sort(df$x), beta = beta,lambda=lambda, delta =delta, mu = mu)
head(df2)

plotar <- df2 %>% select(index,theorical,empirical) %>% melt(id = "index")
tail(plotar)


g3 <- ggplot(df2) +
  geom_line(aes(index,theorical,col='theorical')) + 
  geom_path(aes(index,empirical,col='empirical')) +
  theme(legend.position = c(.20, .85),
        legend.key = element_blank(),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_blank()) +
  scale_colour_manual(values=c(theorical="blue4",empirical="firebrick3")) +
  labs(x="y",y="ECDF")  + guides(color=guide_legend(override.aes=list(fill=NA)))

g3


# tempo de retorno empirico
TT = log(1/(1 - df2$empirical))*100
plot(TT[1:241])


# tempo de retorno teorico
TT = log(1/(1 - df2$theorical))*100 #return period
plot(TT[1:235],sort(fquantil(df2$sobj.surv,beta = beta,lambda=lambda, delta =delta, mu = mu))[1:235]) #return level


head(df2)

df2$quantile <- sort(fquantilm_part(1-df2$sobj.surv,beta = beta,lambda=lambda, delta =delta, mu = mu))

g4 <- ggplot(df2[-nrow(df2),]) +
  geom_point(aes(index,quantile)) + 
  geom_abline(col='red') +labs(x="y")
g4


require(gridExtra)
grid.arrange(g1, g2,g3,g4, ncol=2)

ggsave("imagensR/grid_yellowknife.png", arrangeGrob(g2,g1,g3,g4),
       width = 5, height = 3, dpi = 300, units = "in", device='png')



# tempo de retorno n ----
xp <- seq(.30,.9999,0.0001)
TR <- 1/(1-xp)
zpp <- sort(fquantilm_part(xp,beta = beta,lambda=lambda, delta =delta, mu = mu))
# teorico
zpporiginal <- 100*log(zpp)

# empirico
zppemp <- sort(fquantilm_part(1-df2$sobj.surv,beta = beta,lambda=lambda, delta =delta, mu = mu))
zppemporiginal <- 100*log(zppemp)
xpemp <- 1-df2$sobj.surv
TRemp <- 1/(1-xpemp)

# x temperatura e y em dias, colocar legenda emp e theo, plot original e sep
plot(zpp,TR,type="l",ylim = c(0,2000),xlab="Return level X ",ylab="Return time (days)",main="Yellow Knife")
points(zppemp,pch=16,TRemp,col="red")
legend(1.15, 1600, legend=c("Empirical", "Theorical"),pch = c(19, NA), lty = c(NA, 1),
       col=c("red", "black"), cex=0.8)

# xlim=c(20,31),
plot(zpporiginal,TR,type="l",ylim = c(0,2000),xlab="Return level (C°)",ylab="Return time (days)",main="Yellow Knife")
points(zppemporiginal,pch=16,TRemp,col="red")
legend(20, 1600, legend=c("Empirical", "Theorical"),pch = c(19, NA), lty = c(NA, 1),
       col=c("red", "black"), cex=0.8)



which(xp==.375) #quantil # 2001
TR[751] # tempo de retorno
zpporiginal[751] # temperatura

which(xp==.50) #quantil # 2001
TR[2001] # tempo de retorno
zpporiginal[2001] # temperatura

tail(xp,1000)
which(xp==.9) #quantil # 2001
TR[6001] # tempo de retorno
zpporiginal[6001] # temperatura


which(xp==.95) #quantil
TR[which(xp==.95)] # tempo de retorno
zpporiginal[6501] # temperatura

which(xp==.99) #quantil
TR[which(xp==.99)] # tempo de retorno
zpporiginal[6901] # temperatura

which(xp==.999) #quantil
TR[which(xp==.6991)] # tempo de retorno
zpporiginal[6991] # temperatura


tail(xp,100)

min(zpporiginal)

# SVALBARD ----

teste <- read.csv("dados\\svalbard-climate-1912-2017.csv",encoding="UTF-8",header=T)#,skip=9)
#

library(reshape2)
teste <- melt(teste,id = "YEAR")
svalbar <- teste[teste$value!=999.9,]
svalbar <- svalbar[order(svalbar$YEAR),]

hist(svalbar$value)
x <- (svalbar$value)
# x <- x + abs(min(x))
head(svalbar)
svalbar$x <- exp(x/10)
# transformação pra deixar positiva
x <- exp(x/10)

acf((x))
hist((x))

N<-length(x)
n<-3
tau<-floor(N/n)
M<-matrix(0,tau,1)
j<-1
for (i in 1:tau){
  M[i]<-mean(x[j:(j+n-1)])
  j<-j+n}
head(M)
sum(is.na(M))
hist(M,n=3, prob=T)
lines(density(M),lwd=2)

acf(M)

M<- M[1:35]
length(M)







max(svalbar$value)

beta=  2.09820224047842
lambda=  1.02254548476456
delta =   0.541382456353869
mu=  -1.10301957817666

hist(x, prob = T)
curve(fd_WBm_part(x, beta = beta,lambda=lambda, delta =delta, mu = mu),
      col='green',lwd=2,add=T) 


# LOGVERO -942.7653
# Se dividir por 100 vira: 2199.828
#AIC: AIC = - 2*log L + k * GL, k=2

-2*(-942.7653) + 2*4

#prof vila
chute<-c(1.5,2,0.7)
nlmmodel <- nlm(vero,chute,x=x)
theta_est <- nlmmodel$estimate
-2*((-1)*nlmmodel$minimum) + 2*4


df <- data.frame(x = x,y = fd_WBm_part(x, beta = beta,lambda=lambda, delta =delta, mu = mu),
                 prof=weibivar(x,theta_est),
                 date=c(1:length(x)))

head(df)

g1 <- ggplot(df, aes(x)) + geom_histogram(aes(y = ..density..),fill="gray75", 
                                          binwidth = .15, 
                                          col="black", 
                                          size=.1) +
  geom_line(aes(x,y),col="blue4") +labs(x="y") 

g1



head(df)
# serie temporal de x




g2 <- ggplot(svalbar %>% group_by(YEAR) %>% summarise(x=min(x))) +
  geom_line(aes(YEAR,x),col="blue4") +labs(y="y",x="date")

g2


# Acumuladas
sobj <- survfit(Surv(x,rep(1,length(x)))~1)

dsurv <- data.frame(sobj$time,empirical = 1 - sobj$surv,sobj.surv = sobj$surv)

df2 <- merge(df,dsurv,by.x = "x",by.y = "sobj.time")

df2$index <- sort(df2$x)

df2$theorical <- FD_WBm(sort(df$x), beta = beta,lambda=lambda, delta =delta, mu = mu)

plotar <- df2 %>% select(index,theorical,empirical) %>% melt(id = "index")
tail(plotar)


g3 <- ggplot(df2) +
  geom_line(aes(index,theorical,col='theorical')) + 
  geom_path(aes(index,empirical,col='empirical')) +
  theme(legend.position = c(.20, .85),
        legend.key = element_blank(),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_blank()) +
  scale_colour_manual(values=c(theorical="blue4",empirical="firebrick3")) +
  labs(x="y",y="ECDF")  + guides(color=guide_legend(override.aes=list(fill=NA)))

g3


# tempo de retorno 
TT = 1/df$sobj.surv

head(df2)

df2$quantile <- sort(fquantilm_part(1-df2$sobj.surv,beta = beta,lambda=lambda, delta =delta, mu = mu))

g4 <- ggplot(df2[-nrow(df2),]) +
  geom_point(aes(index,quantile)) + 
  geom_abline(col='red') +labs(x="y")
g4

require(gridExtra)
grid.arrange(g1, g2,g3,g4, ncol=2)

ggsave("imagensR/grid_svalbard.png", arrangeGrob(g2,g1,g3,g4),
       width = 5, height = 3, dpi = 300, units = "in", device='png')



# tempo de retorno n ----
xp <- seq(.3,.9999,0.0001)
TR <- 1/(1-xp)
zpp <- sort(fquantilm_part(xp,beta = beta,lambda=lambda, delta =delta, mu = mu))
# teorico
zpporiginal <- 10*log(zpp)

# empirico
zppemp <- sort(fquantilm_part(1-df2$sobj.surv,beta = beta,lambda=lambda, delta =delta, mu = mu))
zppemporiginal <- 10*log(zppemp)
xpemp <- 1-df2$sobj.surv
TRemp <- 1/(1-xpemp)


length(xp)
plot(zpporiginal)

which(xp==.375) #quantil 
TR[751] # tempo de retorno
zpporiginal[751] # temperatura

which(xp==.50) #quantil # 2001
TR[2001] # tempo de retorno
zpporiginal[2001] # temperatura

xp[6000:6010]
which(xp==.9) #quantil # 2001
TR[6001] # tempo de retorno
zpporiginal[6001] # temperatura

which(xp==.95) #quantil
TR[which(xp==.95)] # tempo de retorno
zpporiginal[6501] # temperatura

which(xp==.99) #quantil
TR[which(xp==.99)] # tempo de retorno
zpporiginal[6901] # temperatura

xp[6990:7000]
which(xp==.9990) #quantil
TR[which(xp==.9991)-1] # tempo de retorno
zpporiginal[6991] # temperatura

which(xp==.9994) #quantil
TR[which(xp==.9994)+1] # tempo de retorno
zpporiginal[2995] # temperatura



# x temperatura e y em dias, colocar legenda emp e theo, plot original e sep
plot(zpp,TR,type="l",ylim = c(0,2000),xlab="Return level X ",ylab="Return time (days)",main="Svalbard")
points(zppemp,pch=16,TRemp,col="red")
legend(1.16, 1600, legend=c("Empirical", "Theorical"),pch = c(19, NA), lty = c(NA, 1),
       col=c("red", "black"), cex=0.8)


plot(zpporiginal,TR,type="l",ylim = c(0,2000),xlab="Return level (C°)",ylab="Return time (days)",main="Svalbard")
points(zppemporiginal,pch=16,TRemp,col="red")
legend(20, 1600, legend=c("Empirical", "Theorical"),pch = c(19, NA), lty = c(NA, 1),
       col=c("red", "black"), cex=0.8)





# DELHI ----


novadelhi <- read.csv("dados\\DailyDelhiClimateTest.csv",encoding="UTF-8",header=T)#,skip=9)
#

x <- novadelhi$meantemp
x <- exp(x/100) 
x <- x[!is.na(x)]

beta=  201.897361761866
lambda=  1.68789121065167
delta =   1.01582450161826
mu=  -1.68984633468315

hist(x, prob = T)
curve(fd_WBm_part(x, beta = beta,lambda=lambda, delta =delta, mu = mu),
      col='green',lwd=2,add=T) 

acf(diff(x))
hist(diff(x))
# LOGVERO -942.7653
#AIC: AIC = - 2*log L + k * GL, k=2

-2*(144.9029) + 2*4

#prof vila
chute<-c(1,1.5,2)
nlmmodel <- nlm(vero,chute,x=x)
nlmmodel
-2*((-1)*nlmmodel$minimum) + 2*4


theta_est <- c(nlmmodel$estimate)#nlmmodel$estimate

df <- data.frame(x = x,y = fd_WBm_part(x, beta = beta,lambda=lambda, delta =delta, mu = mu),
                 prof=weibivar(x,theta_est),
                 date=c(1:length(x)))

head(df)

g1 <- ggplot(df, aes(x)) + geom_histogram(aes(y = ..density..),fill="gray75", 
                                          binwidth = .03, 
                                          col="black", 
                                          size=.1) +
  geom_line(aes(x,y),col="blue4") +labs(x="y") 
g1



# serie temporal de x
g2 <- ggplot(df) +
  geom_line(aes(date,x),col="blue4") + scale_x_continuous(breaks=c(16,45,75,100),
                                                          labels=c("jan","feb","mar","apr")) +labs(y="y")
g2




# Acumuladas
sobj <- survfit(Surv(x,rep(1,length(x)))~1)

dsurv <- data.frame(sobj$time,empirical = 1 - sobj$surv,sobj.surv = sobj$surv)

df2 <- merge(df,dsurv,by.x = "x",by.y = "sobj.time")

df2$index <- sort(df2$x)

df2$theorical <- FD_WBm(sort(df$x), beta = beta,lambda=lambda, delta =delta, mu = mu)

plotar <- df2 %>% select(index,theorical,empirical) %>% melt(id = "index")
tail(plotar)


g3 <- ggplot(df2) +
  geom_line(aes(index,theorical,col='theorical')) + 
  geom_path(aes(index,empirical,col='empirical')) +
  theme(legend.position = c(.20, .85),
        legend.key = element_blank(),
        legend.background=element_rect(fill = alpha("white", 0.3)),
        legend.key.size = unit(0.4, 'cm'),
        legend.title = element_blank()) +
  scale_colour_manual(values=c(theorical="blue4",empirical="firebrick3")) +
  labs(x="y",y="ECDF")  + guides(color=guide_legend(override.aes=list(fill=NA)))

g3

# tempo de retorno 
TT = 1/df$sobj.surv

head(df2)

df2$quantile <- sort(fquantilm_part(1-df2$sobj.surv,beta = beta,lambda=lambda, delta =delta, mu = mu))

g4 <- ggplot(df2[-nrow(df2),]) +
  geom_point(aes(index,quantile)) + 
  geom_abline(col='red') +labs(x="y")
g4


require(gridExtra)
grid.arrange(g1, g2,g3,g4, ncol=2)

ggsave("imagensR/grid_delhi.png", arrangeGrob(g2,g1,g3,g4),
       width = 5, height = 3, dpi = 300, units = "in", device='png')



# tempo de retorno n ----
xp <- seq(.3,.9999,0.0001)
TR <- 1/(1-xp)
zpp <- sort(fquantilm_part(xp,beta = beta,lambda=lambda, delta =delta, mu = mu))
# teorico
zpporiginal <- 100*log(zpp)

# empirico
zppemp <- sort(fquantilm_part(1-df2$sobj.surv,beta = beta,lambda=lambda, delta =delta, mu = mu))
zppemporiginal <- 100*log(zppemp)
xpemp <- 1-df2$sobj.surv
TRemp <- 1/(1-xpemp)

which(xp==.375) #quantil # 2001
TR[751] # tempo de retorno
zpporiginal[751] # temperatura

which(xp==.50) #quantil # 2001
TR[2001] # tempo de retorno
zpporiginal[2001] # temperatura


which(xp==.9000) #quantil # 2001
TR[2001] # tempo de retorno
zpporiginal[2001] # temperatura


which(xp==.95) #quantil
TR[which(xp==.95)] # tempo de retorno
zpporiginal[2501] # temperatura

which(xp==.99) #quantil
TR[which(xp==.99)] # tempo de retorno
zpporiginal[2901] # temperatura

which(xp==.9991) #quantil
TR[which(xp==.9991)-1] # tempo de retorno
zpporiginal[2991] # temperatura

which(xp==.9994) #quantil
TR[which(xp==.9994)+1] # tempo de retorno
zpporiginal[2995] # temperatura



plot(zpp,TR,type="l",ylim = c(0,300),xlab="Return level X ",ylab="Return time (days)",main="Delhi")
points(zppemp,pch=16,TRemp,col="red")
legend(1.28, 250, legend=c("Empirical", "Theorical"),pch = c(19, NA), lty = c(NA, 1),
       col=c("red", "black"), cex=0.8)


plot(zpporiginal,TR,type="l",ylim = c(0,300),xlab="Return level (C°)",ylab="Return time (days)",main="Delhi")
points(zppemporiginal,pch=16,TRemp,col="red")
legend(24.5, 250, legend=c("Empirical", "Theorical"),pch = c(19, NA), lty = c(NA, 1),
       col=c("red", "black"), cex=0.8)








# MODELO DE REGRESSÃO YELLOWKNIFE ----

arquivos <- list.files(path = "D:\\foque\\unb\\!Mestrado\\Mono\\R\\dados\\yellowknifeanos",full.names=T)

yknife <- read.csv(arquivos[1],encoding="UTF-8",header=T)
for(i in 2:length(arquivos)){
  a <- read.csv(arquivos[i],encoding="UTF-8",header=T)
  yknife <- rbind(yknife,a)
  
}

head(yknife)
tail(yknife)

summary(yknife$Min.Temp...C.)

yknife$Season <- factor(yknife$Month, levels = c(1,2,3,4,5,6,7,8,9,10,11,12), 
                        labels = c("Winter","Winter","Winter", "Spring", "Spring", "Spring",
                                   "Summer","Summer","Summer","Fall","Fall","Fall"))



x <- yknife$Min.Temp...C.
hist(x)
x <- exp(x/100) # muda ser sobre 10 ou 100
yknife <- yknife[!is.na(x),]
x <- x[!is.na(x)]
seas <- yknife$Season

N<-length(x)
n<-90
tau <- floor(N/n)
M<-matrix(0,tau,1)
Mseason <- c()
j<-1
for (i in 1:tau){
  M[i]<-min(x[j:(j+n-1)])
  Mseason[i] <- seas[j:(j+n-1)][which.min(x[j:(j+n-1)])]
  j<-j+n}
head(M)
sum(is.na(M))
hist(M, prob=T)
lines(density(M),lwd=2)
x <- M

length(x)/4

Mseason <- as.factor(Mseason)
levels(Mseason) <-c("Winter","Spring","Summer","Fall")


# determinando as covariaveis

yknifereg <- data.frame(x,season=Mseason)


yknifereg$x0 <- rep(1,nrow(yknifereg))

yknifereg$x1 <- as.numeric(yknifereg$season=="Winter")
yknifereg$x2 <- as.numeric(yknifereg$season=="Spring")
yknifereg$x3 <- as.numeric(yknifereg$season=="Fall")


x <- yknifereg$x
x0 <- yknifereg$x0 
x1 <- yknifereg$x1 
x2 <- yknifereg$x2 
x3 <- yknifereg$x3 


# log-vero regressão 
lvero_reg_1 <- function(par){
  beta0 <- par[1]
  beta1 <- par[2]
  beta2 <- par[3]
  beta3 <- par[4]
  
  lambda <- exp(-(beta0*x0+beta1*x1+beta2*x2+beta3*x3))

  beta <-par[5]
  delta<-par[6]
  mu <-par[7]
  m <- (-mu)^(1/(delta+1))
  v = x - m
  
  # L1
  ind1 <- ifelse(x<m,1,0)

  #L2
  ind2 <- ifelse(x>m,1,0)
  
  
  # L1
  options(warn=-1) 
  a1 <- log(beta/lambda) +log(delta+1)
  a2 <- - ((-(-v)^(delta+1)-mu)/lambda)^beta
  a3 <- (beta-1)*(log((-(-v)^(delta+1)-mu)/lambda))
  a4 <- delta*log(-v)
  a = a1+a2+a3+a4
  
  #L2
  b1 <- log(beta/lambda) +log(delta+1) 
  b2 <- - ((v^(delta+1)-mu)/lambda)^beta
  b3 <- (beta-1)*(log((v^(delta+1)-mu)/lambda))
  b4 <- delta*log(v)
  b = b1+b2+b3+b4
  
  # Coloca 0 no NaN (talvez nem precisasse disso,
  # só remover e somar)
  b[is.nan(b)] <- 0   
  a[is.nan(a)] <- 0  
  b[(b==-Inf)] <- 0   
  a[a==-Inf] <- 0  
  
  
  f = sum(ind1*a,ind2*b)
  
  
  # print(a)
  # print(b)
  # print(ind1)
  # print(ind2)
  return(f)#sum(f))
}

df <- cbind(runif(1000,0,1),runif(1000,0,1),runif(1000,0,1),runif(1000,0,1),
            runif(1000,0,30),runif(1000,0.03,10),-runif(1000,0.5,5))

valormax <- c()
estimados <- list()
vari <- list()
codigo <- c()
confirm <- c()

# estimação
for(i in 1:nrow(df)){
  tryCatch({
    a <- optim(df[i,] , lvero_reg_1,
               control=list(fnscale=-1,maxit=10000),
               lower=c(rep(-Inf,2),0.0001,0.06,-Inf), 
               upper=c(rep(Inf,6),-0.01),
               method="L-BFGS-B",
               hessian=T) 
    valormax[i] <- a$value
    estimados[[i]] <- a$par
    vari[[i]] <- solve((-1)*a$hessian)
    confirm[i] <- any(diag(solve((-1)*a$hessian))<0)
    codigo[i] <- a$convergence
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  cat(i,a$value,"\n")
}

# valores max
head(sort(valormax,decreasing=T),10)


# valores dos parametros 
dest <- head(estimados[order(valormax,decreasing = T)],20)
quemsao <- head(df[order(valormax,decreasing = T),],20)

# confirmação de convergencia e hessiana
head(codigo[order(valormax,decreasing = T)],20)
head(confirm[order(valormax,decreasing = T)],20)


# pegando o valor dos parâmetros estimados
dest <- dest[[1]]


# winter 
lambda1 <- exp(-(dest[1]+dest[2]))
beta1 = dest[5]
delta1 =dest[6]
mu1 = dest[7]

xwin <- yknifereg[yknifereg$season=="Winter",]$x

hist(xwin, prob = T)
points(seq(0.5,2,0.01),fd_WBm_part(seq(0.5,2,0.01), beta = beta1,lambda=lambda1, delta = delta1, mu = mu1),
       col='green',type='l',lwd=2)


df <- data.frame(x = xwin,y = fd_WBm_part(xwin, beta = beta1,lambda=lambda1, delta =delta1, mu = mu1),
                 date=c(1:length(xwin)))

library(ggplot2)
theme_set(theme_bw())
g1 <- ggplot(df, aes(x)) + geom_histogram(aes(y = ..density..),fill="gray75",
                                          bins = 5,
                                          #binwidth = .009,
                                          col="black",
                                          size=.1) +
  labs(x="x:winter") +
  geom_line(aes(x,y),col="black")

g1
ggsave(filename = paste("imagensR/regwinter.png",sep=""),g1,
       width = 5, height = 3, dpi = 300, units = "in", device='png')



# spring
lambda2 <- exp(-(dest[1]+dest[3]))
beta2 = dest[5]
delta2 =dest[6]
mu2 = dest[7]


xspr <- yknifereg[yknifereg$season=="Spring",]$x #- .2

hist(xspr, prob = T)
points(seq(0.6,1.05,0.01),fd_WBm_part(seq(0.6,1.05,0.01), beta = beta2,lambda=lambda2, delta = delta2, mu = mu2),
       col='green',type='l',lwd=2)


df <- data.frame(x = xspr,y = fd_WBm_part(xspr, beta = beta2,lambda=lambda2, delta =delta2, mu = mu2),
                 date=c(1:length(xspr)))
theme_set(theme_bw())
g2 <- ggplot(df, aes(x)) + geom_histogram(aes(y = ..density..),fill="gray75",
                                          bins=9,
                                          #binwidth = .05,
                                          col="black",
                                          size=.1) +
  labs(x="x:spring") +
  geom_line(aes(x,y),col="black")

g2
ggsave(filename = paste("imagensR/regspring.png",sep=""),g2,
       width = 5, height = 3, dpi = 300, units = "in", device='png')


# Autumn

lambda3 <- exp(-(dest[1]+dest[4]))
beta3 = dest[5]
delta3 =dest[6]
mu3 = dest[7]

xfal <- yknifereg[yknifereg$season=="Fall",]$x

hist(xfal, prob = T)
points(seq(0.6,1.2,0.01),fd_WBm_part(seq(0.6,1.2,0.01), beta = beta3,lambda=lambda3, delta = delta3, mu = mu3),
       col='green',type='l',lwd=2)


xfal2 <- append(xfal,seq(max(xfal),1.25,0.01))


df <- data.frame(x = xfal,y = fd_WBm_part(xfal, beta = beta3,lambda=lambda3,
                                           delta =delta3, mu = mu3),
                 date=c(1:length(xfal)))
theme_set(theme_bw())
g3 <- ggplot(df, aes(x)) + geom_histogram(aes(y = ..density..),fill="gray75",
                                          #binwidth = .05,
                                          bins=5,
                                          col="black",
                                          size=.1) + #xlim(c(0.6,1)) +
  labs(x="x:autumn") +
  geom_line(aes(x,y),col="black")

g3
ggsave(filename = paste("imagensR/regfall.png",sep=""),g3,
       width = 5, height = 3, dpi = 300, units = "in", device='png')



# summer

lambda4 <- exp(-(dest[1]))
beta4 = dest[5]
delta4 =dest[6]
mu4 = dest[7]

xsum <- yknifereg[yknifereg$season=="Summer",]$x


hist(xsum, prob = T)
points(seq(0.7,2,0.01),fd_WBm_part(seq(0.7,2,0.01), beta = beta4,lambda=lambda4, delta = delta4, mu = mu4),
       col='green',type='l',lwd=2)


df <- data.frame(x = xsum,y = fd_WBm_part(xsum, beta = beta4,lambda=lambda4, delta =delta4, mu = mu4),
                 date=c(1:length(xsum)))
theme_set(theme_bw())
g4 <- ggplot(df, aes(x)) + geom_histogram(aes(y = ..density..),fill="gray75", 
                                          #binwidth = .05, 
                                          bins = 8,
                                          col="black", 
                                          size=.1) +
  labs(x="x:summer") +
  geom_line(aes(x,y),col="black")

g4
ggsave(filename = paste("imagensR/regsummer.png",sep=""),g4,
       width = 5, height = 3, dpi = 300, units = "in", device='png')



ggsave("imagensR/regall0.png", arrangeGrob(g1, g2,g3,g4),
       width = 5, height = 3, dpi = 300, units = "in", device='png')

dest
