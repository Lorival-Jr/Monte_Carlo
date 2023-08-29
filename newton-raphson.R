U         <- function(x, par)
{    # Vetor escore é um vetor que contém a primeira derivada em relação a cada parâmetro, nesse caso é a derivada primeira em relação a theta, derivada primeira em relação a rho
  theta <- par[1]
  rho   <- par[2]
  n       <- length(x)
  
  dtheta  <- 		  -sum(((2*(x*((x*theta)/(theta+1)+1)*exp(-(x*theta))*rho-(x/(theta+1)-(x*theta)/(theta+1)^2)*exp(-(x*theta))*rho))/(1-((x*theta)/(theta+1)+1)*exp(-(x*theta))*rho)))+n*(2/theta-1/(theta+1))-sum(x)
  

  
  drho1   <- ((x*theta)/(theta+1)+1)
  drho2   <- exp(-(x*theta))
  drho    <- sum((2*drho1*drho2)/(1-drho1*drho2*rho))-n/(1-rho)
  
  vetor   <- matrix(nrow = 2, c(dtheta, drho))
  
  return(vetor)
}

set.seed(123)
U(x, c(2,0.3))




dtheta2   <- function(x, par)           # Derivada segunda da log-verossimilhança em relação a theta
{
  n       <- length(x)
  theta   <- par[1]
  rho     <- par[2]
  
  e      <- exp(-x*theta)
  p1     <- x/(par[1]+1)
  p2     <- p1*theta
  p3     <- p2/(theta + 1)
  
  d1      <- sum((2*(x*(p2+1)*e*par[2]-(p1-p3)*e*par[2])^2)/(1-(p2+1)*e*par[2])^2)
  d2      <- -sum((2*(-(x^2*(p2+1)*e*par[2]))+2*x*(p1-p3)*e*par[2]))
  d3      <- -sum(((2*x*par[1])/(par[1]+1)^3-(2*x)/(par[1]+1)^2*e*par[2])/(1-(p2+1)*e*par[2]))+n/(par[1]+1)^2-(2*n)/par[1]^2
  dd      <-  d1 + d2 + d3
  
  return(dd)
}
set.seed(123)
dtheta2(rlingley_geom(1000, c(3,0.5)), par = c(3, 0.5))

# Rho-Rho

drho2     <- function(x, par)          # Derivada segunda da log-verossimilhança em relação a rho
{
  n       <- length(x)
  theta   <- par[1]
  rho     <- par[2]
  
  
  p1      <- ((x*theta)/(theta+1)+1)
  
  dd      <- sum((2*p1^2*exp(-(2*x*par[1])))/(1-p1*exp(-(x*par[1]))*par[2])^2)-n/(1-par[2])^2
  
  return(dd) 
}
set.seed(123)
drho2(rlingley_geom(1000, c(3,0.5)), par = c(3, 0.5))

# Theta-Rho

dthetarho <- function(x, par)          # Derivada segunda da log-verossimilhança em relação a theta e depois a rho
{
  theta   <- par[1]
  rho     <- par[2]
  
  e  <- exp(-(x*par[1]))
  
  p0 <- (x)/(theta + 1)
  p1 <- p0*theta + 1
  p2 <- p0*theta/(theta + 1)
  d1 <- -sum(((2*(x*p1*e-(x/(par[1]+1)-p2)*e))/(1-p1*e*par[2])))
  d2 <- -sum((2*p1*e*(x*p1*e*par[2]-(x/(par[1]+1)-p2)*e*par[2]))/(1-p1*e*par[2])^2)
  dd <-  d1 + d2
  


  
  return(dd)
}
set.seed(123)
dthetarho(rlingley_geom(1000, c(3,0.5)), par = c(3, 0.5))



J <- function(x, par)
{
  
  jacobiano      <- matrix(nrow = 2, ncol = 2)
  jacobiano[1,1] <- dtheta2(x, par)
  jacobiano[1,2] <- dthetarho(x, par)
  jacobiano[2,1] <- jacobiano[1,2]
  jacobiano[2,2] <- drho2(x, par)
  
  return(jacobiano)
  
}

J <- function(x, par)
{
  theta   <- par[1]
  rho     <- par[2]
  n       <- length(x)
  
  e  <- exp(-(x*theta))
  
  p0 <- (x)/(theta + 1)
  p1 <- p0*theta + 1
  p2 <- p0*theta/(theta + 1)
  d1 <- -sum(((2*(x*p1*e-(x/(par[1]+1)-p2)*e))/(1-p1*e*par[2])))
  d2 <- -sum((2*p1*e*(x*p1*e*par[2]-(x/(par[1]+1)-p2)*e*par[2]))/(1-p1*e*par[2])^2)
  
  dthetarho <- d1 + d2 
  ###################
  
  drho2      <- sum((2*p1^2*exp(-(2*x*par[1])))/(1-p1*exp(-(x*par[1]))*par[2])^2)-n/(1-par[2])^2
  ###################
  
  d1      <- sum((2*(x*(p1)*e*par[2]-(p0-p2)*e*par[2])^2)/(1-(p1)*e*par[2])^2)
  d2      <- -sum((2*(-(x^2*(p1)*e*par[2]))+2*x*(p0-p2)*e*par[2]))
  d3      <- -sum(((2*x*par[1])/(par[1]+1)^3-(2*x)/(par[1]+1)^2*e*par[2])/(1-(p1)*e*par[2]))+n/(par[1]+1)^2-(2*n)/par[1]^2
  dtheta2 <-  d1 + d2 + d3
  
  jacobiano      <- matrix(nrow = 2, ncol = 2)
  jacobiano[1,1] <- dtheta2
  jacobiano[1,2] <- dthetarho
  jacobiano[2,1] <- jacobiano[1,2]
  jacobiano[2,2] <- drho2
  
  return(jacobiano)
  
}
set.seed(9999)
J(rlingley_geom(1000, c(2,0.5)), c(2,0.5))
J(rlingley_geom(1000, c(0.3, 0.5)), c(0.3, 0.5))





par <- matrix(c(6, 0.5), nrow = 2)                # Os valores que o Newton-raphson precisa
set.seed(9999)
x <-  rlingley_geom(50000, c(6, 0.5))
itr <- 5000
erro <- 0.0000001
i <- 1

while(i < itr)                                    # Método de Newton-Raphson, mas tá explodindo
{                                                 # Erro pode ser nas derivas, na verossimilhança 
  # ou na própria implementação do método
  i <- i + 1
  
  par_novo <- par - solve(J(x, par)) %*% U(x, par)
  
  if(max(abs(U(x, par_novo) - U(x, par))) < erro) break # Critério de parada
  
  par <- par_novo
  
  cat('Iteração:', i, 'theta:', par_novo[1], 'rho:', par_novo[2], "\n")
  
  
  
}



itr <- 500
erro <- 10e-6
i <- 1


par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000

simulacoes_newton <- array(c(rep(0,6)), dim=c(N,15,9,10))
simulacoes_newton

set.seed(9999)
set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{ 
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9)                         # Combinação de parâmetros
  { 
  par0 <- par_comb[[index_par]]
  amostra <- rlingley_geom(n, par=par0)               # Amostra
 
  iteracao <- 1
  
  par <- par0 
  while(iteracao < itr)                                    # Método de Newton-Raphson, mas tá explodindo
  {                                                 # Erro pode ser nas derivas, na verossimilhança 
    # ou na própria implementação do método
    iteracao <- iteracao + 1
    
    variancia <- try(- solve(J(amostra, par)),F)
    if(typeof(variancia) == 'character') {variancia <- c(NA, NA); break}
    
    par_novo <- par + variancia %*% U(amostra, par)
    
    if(is.nan(U(amostra, par_novo))){break}
    
    if(max(abs(U(amostra, par_novo) - U(amostra, par))) < erro) break # Critério de parada
    
    par <- par_novo
    
    #cat('Iteração:', i, 'theta:', par_novo[1], 'rho:', par_novo[2], "\n")
    
    
    
  }
  
  valores <- c(par_novo[1], par_novo[2], par0[1], par0[2], n, variancia[1], variancia[4], rep(0,8))
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
  
  cat('itr:', i, '-' , valores, '\n')
  simulacoes_newton[i, ,index_par, index_n] <- valores
  
  }
  }
  
}

save(simulacoes_newton, file='simulacoes_newton.RData')



dlindley_geom <- function(x, par)
{
  e           <- exp(-par[1] * x) 
  p1          <- (par[1]^2 / (par[1] + 1)) * (1-par[2]) * (1+ x) * e
  p2          <- (1 - par[2]*(1 + (par[1] * x)/(par[1] + 1)) * e)^(-2)
  pdf         <- p1 * p2
  
  return(pdf)
}
# restrições theta > 0 e 0 < p < 1
#  stopifnot(p < 1, p > 0, theta > 0, x > 0) no R

plindley_geom <- function(x, par) 
{
  e              <- exp(-par[1] * x)
  p1             <- (1 + (par[1]*x)/(par[1] + 1))*e
  cdf            <- (1 - p1)/(1 - par[2]*p1) 
  
  return(cdf)
}

#install.packages('LambertW')
library('LambertW')
??LambertW

qlindley_geom <- function(u, par)
  
{
  e           <- exp(-par[1] - 1)
  p1          <- ((u - 1) * (par[1] + 1) * e) / (1 - par[2]*u)
  quantile    <- -1 - 1 / par[1] - W(p1, -1)* par[1]^(-1)
  
  return(quantile)
}

rlingley_geom <- function(n, par)
{
  
  u <- runif(n, min = 0, max = 1) 
  rnd.values <- qlindley_geom(u, par) 
  return(rnd.values)
  
}












#install.packages('parallel')
#install.packages("doParallel")
library(plyr)




nr <- function(i, itr = 5000, erro = 10e-6, N = 50000, par_comb = list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), 
                                                                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                                                                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))){
  

    for (index_n in 1:10)                         # Tamanho da amostra
    { n <- seq(10, 100, 10)[index_n]
    
    for (index_par in 1:9)                         # Combinação de parâmetros
    { 
      par0 <- par_comb[[index_par]]
      amostra <- rlindley_geom(n, par=par0)               # Amostra
      
      iteracao <- 1
      
      par <- par0 
      while(iteracao < itr)                                    # Método de Newton-Raphson, mas tá explodindo
      {                                                 # Erro pode ser nas derivas, na verossimilhança 
        # ou na própria implementação do método
        iteracao <- iteracao + 1
        
        variancia <- try(- solve(J(amostra, par)),F)
        if(typeof(variancia) == 'character') {variancia <- c(NA, NA); break}
        
        par_novo <- par + variancia %*% U(amostra, par)
        
        if(is.nan(U(amostra, par_novo))){break}
        
        if(max(abs(U(amostra, par_novo) - U(amostra, par))) < erro) break # Critério de parada
        
        par <- par_novo
        
        #cat('Iteração:', i, 'theta:', par_novo[1], 'rho:', par_novo[2], "\n")
        
        
        
      }
      
      valores <- c(par_novo[1], par_novo[2], par0[1], par0[2], n, variancia[1], variancia[4], rep(0,8))
      # Valores recebe o que queremos dessa bagaça toda,
      # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
      
      cat('itr:', i, '-' , valores, '\n')
      simulacoes_newton[i, ,index_par, index_n] <<- valores
      
    
    }
    
   }
}

nr(1:100)

library(plyr)
library(doParallel)
cl <- makeCluster(2)
registerDoParallel(cl)
opts <- list(preschedule=TRUE)
clusterSetRNGStream(cl, 123)
clusterCall(cl, function() library(LambertW))
clusterExport(cl, c("rlindley_geom", 'qlindley_geom', 'J', 'U', 'simulacoes_newton'))

simulacoes_newton <- array(c(rep(0,6)), dim=c(N,15,9,10))
simulacoes_newton

library(LambertW)
llply(1:50000, nr, .parallel = T, .paropts = list(.options.snow=opts), )

help(llply)
parallel::parLapply(cl,
                    1:50000,
                    nr)

parallel::stopCluster(cl)
