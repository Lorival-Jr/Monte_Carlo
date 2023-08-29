# Geração de Valores Aleatórios -------------------------------------------


# 9. Geração de Valores Aleatórios

# Lindley Geométrica ------------------------------------------------------

### Forma analítica

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

rlindley_geom <- function(n, par)
{
  
  u <- runif(n, min = 0, max = 1) 
  rnd.values <- qlindley_geom(u, par) 
  return(rnd.values)
  
}
## Testando

#install.packages("LindleyPowerSeries")


# Pagina 4
# https://cran.r-project.org/web/packages/LindleyPowerSeries/LindleyPowerSeries.pdf
help(dlindleygeometric)

dlindley_geom(c(1, 1.5, 4), c(1, 0.4))

plindley_geom(c(1, 1.5, 4), c(1, 0.4))

qlindley_geom(c(0.3,0.5,0.8), c(2, 0.4))

rlingley_geom(100, c(2, 0.4))

hist(rlingley_geom(1000, c(0.8, 0.4)), freq = F)

par(mfrow = c(1,2))
hist(rlingley_geom(1000, c(0.2, 0.4)), freq = F)
hist(rgeom(1000, 0.4), freq = F)
par(mfrow = c(1,1))


### Método numérico

dlindley_geom <- function(x, theta, p)
{
  e           <- exp(-theta * x) 
  p1          <- (theta^2 / (theta + 1)) * (1-p) * (1+ x) * e
  p2          <- (1 - p*(1 + (theta * x)/(theta + 1)) * e)^(-2)
  pdf         <- p1 * p2
  
  return(pdf)
}
# restrições theta > 0 e 0 < p < 1

plindley_geom_num   <-      function(x, theta, p)
{
  my.int            <-      function(x, theta, p) # Criando uma função pra utilizar no apply
  { # Ela integra a fdp, recebendo o parametro lambda pra fdp e o lower sendo o intervalo de integração inf = 0 e o sup = x
    integrate(dlindley_geom, theta = theta, p = p, lower = 0, upper = x)$value # é pego do resultado apenas o valor
  } 
  
  sapply(x, FUN=my.int, theta = theta, p = p) # Aplica nos dados a integral da fdp, ou seja a fda
}


qlindley_geom_num  <-      function(q, theta, p, lower=0, upper) # Função quantil
{ 
  f                 <-      function(P, fixed)  # Função pra ser aplicada no unirrot                    
  {
    lambda          <-      fixed$theta        # fixed é uma lista com 3 dimensões
    p               <-      fixed$p
    q               <-      fixed$q             
    criterion       <-      q - plindley_geom(P, theta, p)   # É o inverso da probabiladade da FDA                     
    return(criterion)
  }
  P                 <-      numeric(length(q)) # Para cada q passado vai achar uma raiz                                     
  for(i in 1:length(q))
  {
    fixed           <-      list(theta = theta, p = p, q = q[i]) # É a lista que será passado para função f
    root.p          <-      uniroot(f, lower = lower, upper = upper, fixed = fixed) 
    # unirrot procura as raizes da função f, com o intervalo de integração de lower a upper, o fixed é a lista que a f recebe
    P[i]            <-      root.p$root # pega apenas a raiz pra retornar
  }
  return(P)
}

## Testando 
plindley_geom_num(c(1, 1.5, 4), 1, 0.5)
plindley_geom(c(1, 1.5, 4), 1, 0.5)

qlindley_geom_num(c(0.3,0.5,0.8), 2, 0.4, 0, 100)
qlindley_geom(c(0.3,0.5,0.8), 2, 0.4)


log_lindley_geometrica <- function(x, par) # par[1] será theta, par[2] é o p
{  #Restrições: par[1] > 0 e 0 < par[2] < 1
  
  
  # 2*n*log(par[1]) + n*log(1-par[2]) + sum(log(1+x)) + 
  #   sum(log(exp(-theta[1]*x))) - n*log(par[1] + 1) -
  #   2*sum(log(1-theta[2]*(1+(theta[1]*x))))
  
  n  <- length(x)
  f1 <- 2*n*log(par[1])
  f2 <- n*log(1-par[2])
  f3 <- sum(log(1+x))
  f4 <- sum(-par[1]*x)
  f5 <- n*log(par[1]+1)
  f6 <- 2*sum(log(1-par[2]*(1+par[1]*x/(par[1]+1))*exp(-x*par[1])))
  
  ll <- f1 + f2 + f3 + f4 - f5 - f6
  return(ll)
  
}

amostra <- rlingley_geom(10000, 5, 0.8)

optim(par = c(1, 0.3), fn = log_lindley_geometrica, x = amostra, control = list(fnscale = -1))
#Restrições: par[1] > 0 e 0 < par[2] < 1


# Escore de Fisher --------------------------------------------------------

# Derivadas segundas

# Essas derivadas são usadas no Jacobiano

# Theta-Theta

dtheta2   <- function(x, par)           # Derivada segunda da log-verossimilhança em relação a theta
{
  n       <- length(x)
  f1      <- (x*par[1])/(par[1] + 1) + 1
  f2      <- exp(-x*par[1])*par[2]
  f3      <- x/(par[1] + 1)
  f4      <- (x*par[1])/((par[1] + 1)^2)
  
  d1      <- 2*(((x*(f1)*f2) - (f3 - f4)*f2)^2)/(1 - (f1)*f2)^2
  d2      <- (2*(-(x^2*(f1)*f2) + 2*x*(f3 - f4)*f2 - ((2*x*par[1])/(par[1] + 1)^3 - 2*x/(par[1] + 1)^2)*f2))/(1 - f1*f2)
  
  dd      <- sum(d1 - d2) + (n/(par[1]+1)^2 - 2*n/(par[1]^2))
  
  return(dd)
}
set.seed(123)
dtheta2(rlingley_geom(1000, c(3,0.5)), par = c(3, 0.5))

# Rho-Rho

drho2     <- function(x, par)          # Derivada segunda da log-verossimilhança em relação a rho
{
  n       <- length(x)
  f1      <- ((x*par[1])/(par[1]+1)+1)
  
  d1      <- 2*f1^2*exp(-2*x*par[1])
  d2      <- (1 - f1*exp(-x*par[1])*par[2])^2
  
  dd      <- sum(d1/d2) -  n/(1 - par[2])^2
  
  return(dd) 
}
set.seed(123)
drho2(rlingley_geom(1000, c(3,0.5)), par = c(3, 0.5))

# Theta-Rho

dthetarho <- function(x, par)          # Derivada segunda da log-verossimilhança em relação a theta e depois a rho
{
  f1      <- ((x*par[1])/(par[1]+1)+1)
  f2      <- exp(-x*par[1])
  f3      <- x/(par[1] + 1)
  f4      <- (x*par[1])/((par[1] + 1)^2)
  f5      <- f2*par[2]
  
  d1 <- (2*x*f1*f2)/(1-f1*f5)
  
  d2 <- (2*(f3-f4)*f2)/(1-f1*f5)
  
  d3 <- ((2*f1*f2*(x*f1*f5-(f3-f4)*f5))/(1-f1*f5)^2)
  
  dd <- sum(-d1 + d2 - d3)
  
  
  return(dd)
}
set.seed(123)
dthetarho(rlingley_geom(1000, c(3,0.5)), par = c(3, 0.5))


### Implementação do vetor escore

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


simulacoes_newton[,8,,]  <-  simulacoes_newton[,1,,] - simulacoes_newton[,3,,]     # Vício de rho
simulacoes_newton[,9,,]  <- (simulacoes_newton[,1,,] - simulacoes_newton[,3,,])^2 # EQM de rho
simulacoes_newton[,10,,] <- (simulacoes_newton[,6,,])^0.5
simulacoes_newton[,11,,] <- ((simulacoes_newton[,1,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_newton[,6,,])) < simulacoes_newton[,3,,]) & ((simulacoes_newton[,1,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_newton[,6,,])) > simulacoes_newton[,3,,]) # Prob cobertura rho

simulacoes_newton[,12,,] <-  simulacoes_newton[,2,,] - simulacoes_newton[,4,,]
simulacoes_newton[,13,,]  <- (simulacoes_newton[,2,,] - simulacoes_newton[,4,,])^2
simulacoes_newton[,14,,] <- (simulacoes_newton[,7,,])^0.5
simulacoes_newton[,15,,] <- ((simulacoes_newton[,2,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_newton[,7,,])) < simulacoes_newton[,4,,]) & ((simulacoes_newton[,2,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_newton[,7,,])) > simulacoes_newton[,4,,]) # Prob cobertura rho

simulacoes_newton

diagnostico_newton <- simulacoes_newton[,c(5:15),,]

save(simulacoes_newton,  file = 'simulacoes_newton.Rdata')
save(diagnostico_newton, file = 'diagnostico_newtonr.Rdata')

### OPTIM - Simulações

# Número de simulações (N): 50000.
# Tamanhos de amostra (n): 10, 20, ..., 100.
# Valores paramétricos: devem ser consideradas combinações a fim de gerar nove cenários.

# Métodos - Newton-Raphson, Nelder-Mead, BFGS, CG, L-BFGS-B, SANN

help(optim)



# Nelder-Mead -------------------------------------------------------------
library(LambertW)


# Uma tentativa de maximização com optim

optim(par = c(1, 0.3), fn = log_lindley_geometrica, x = rlingley_geom(1000, c(4, 0.1)), control = list(fnscale = -1), method = 'Nelder-Mead')
#Restrições: par[1] > 0 e 0 < par[2] < 1


# oq queremos guardar? par, convergencia, par pop, n, e as variâncias que vem da hessiana

par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8),   # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),   # Daí são as combinações
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000                                                # O N de simulações pedido é 50000

simulacoes_nelder <- array(c(rep(0,6)), dim=c(N,15,9,10))  # Esse array vai guardar os resultados
simulacoes_nelder                                         # Basicamente são 90 matrizes dentro de um array
 
dim(simulacoes_nelder) # Serão 50 mil linhas, 15 colunas e 90 matrizes

# As dimensões representam, respectivamente, Linha, coluna, dimensão referente a comb dos parâmetros, dimensão do tamanho de amostra
set.seed(9999)
for (i in 1:N) # Número de simulações
{ 
  for (index_n in 1:10) # Tamanho da amostra
  {n <- seq(10, 100, 10)[index_n]
  
    for (index_par in 1:9) # Combinação de parâmetros
    { par <- par_comb[[index_par]]
      amostra <- rlingley_geom(n, par=par)     # Amostra
      op      <- try(optim(par = par, # Chute inicial
                fn = log_lindley_geometrica,   # Log-Verossimilhança
                x = amostra,                   # Amostra
                control = list(fnscale = -1),
                method = 'Nelder-Mead',        # Método
                hessian = T), T)                  # Calcular a hessiana
      
      if(typeof(op) == 'character')
      { op      <- try(optim(par = par, # Chute inicial
                             fn = log_lindley_geometrica,   # Log-Verossimilhança
                             x = amostra,                   # Amostra
                             control = list(fnscale = -1),
                             method = 'Nelder-Mead',        # Método
                             hessian = F)) 
      
      }
      
      h <- try(solve(op$hessian))              # Tenta inverter a hessiana
      if(typeof(h) == 'character') {h <- c(NA, NA, NA, NA)}  # Se não for invetível, ele guarda o erro em character
                                                             # Daí se o tipo for character, h vira um vetor de NA
      
      valores <- c(op$par[1], op$par[2], par[1], par[2], n, -h[1], -h[4], rep(0,8))
      # Valores recebe o que queremos dessa bagaça toda,
      # theta_estimado, rho_estimado, theta_real, rho_real, n, variância_rho, variância_theta, e os zeros serão substituídos fora do for por vício_theta, eqm_theta, erro padrão_theta, probabilidade de cobertura_theta, vício_rho, eqm_rho, erro padrão_rho, probabilidade de cobertura_rho
      
      cat('itr:', i, '-' , valores, '\n')    # Inútil, é só pra vc ficar vendo oq ta acontecendo
      
      simulacoes_nelder[i, ,index_par, index_n] <- valores  # Guarda na tabela
      
    }
  }
  
}

simulacoes_nelder #theta_estimado, rho_estimado, theta_real, rho_real, n, variância_rho, variância_theta, e os zeros serão substituídos fora do for por vício_theta, eqm_theta, erro padrão_theta, probabilidade de cobertura_theta, vício_rho, eqm_rho, erro padrão_rho, probabilidade de cobertura_rho


# As dimensões representam, respectivamente, Linha, coluna, dimensão referente a comb dos parâmetros, dimensão do tamanho de amostra
simulacoes_nelder[,,9,1]


simulacoes_nelder[,8,,]  <-  simulacoes_nelder[,1,,] - simulacoes_nelder[,3,,]     # Vício de rho
simulacoes_nelder[,9,,]  <- (simulacoes_nelder[,1,,] - simulacoes_nelder[,3,,])^2 # EQM de rho
simulacoes_nelder[,10,,] <- (simulacoes_nelder[,6,,])^0.5
simulacoes_nelder[,11,,] <- ((simulacoes_nelder[,1,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_nelder[,6,,])) < simulacoes_nelder[,3,,]) & ((simulacoes_nelder[,1,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_nelder[,6,,])) > simulacoes_nelder[,3,,]) # Prob cobertura rho

simulacoes_nelder[,12,,] <-  simulacoes_nelder[,2,,] - simulacoes_nelder[,4,,]
simulacoes_nelder[,13,,]  <- (simulacoes_nelder[,2,,] - simulacoes_nelder[,4,,])^2
simulacoes_nelder[,14,,] <- (simulacoes_nelder[,7,,])^0.5
simulacoes_nelder[,15,,] <- ((simulacoes_nelder[,2,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_nelder[,7,,])) < simulacoes_nelder[,4,,]) & ((simulacoes_nelder[,2,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_nelder[,7,,])) > simulacoes_nelder[,4,,]) # Prob cobertura rho

simulacoes_nelder

diagnostico_nelder <- simulacoes_nelder[,c(5:15),,]

save(simulacoes_nelder,  file = 'simulacoes_nelder.Rdata')
save(diagnostico_nelder, file = 'diagnostico_nelder.Rdata')


# BFGS --------------------------------------------------------------------

library(LambertW)

# Uma tentativa de maximização com optim
optim(par = c(1, 0.3), fn = log_lindley_geometrica, x = rlingley_geom(1000, c(4, 0.1)), control = list(fnscale = -1), method = 'BFGS', hessian = T)
 
par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000

simulacoes_bfgs <- array(c(rep(0,6)), dim=c(N,15,9,10))
simulacoes_bfgs

dim(simulacoes_bfgs) # Serão 50 mil linhas, 8 colunas e 90 matrizes
# As dimensões representam, respectivamente, Linha, coluna, dimensão referente a comb dos parâmetros, dimensão do tamanho de amostra


set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
    for (index_par in 1:9)                         # Combinação de parâmetros
    { par <- par_comb[[index_par]]
    amostra <- rlingley_geom(n, par=par)           # Amostra
    op      <- try(optim(par = par,            # Chute inicial
                     fn = log_lindley_geometrica,  # Log-Verossimilhança
                     x = amostra,                  # Amostra
                     control = list(fnscale = -1),
                     method = 'BFGS',              # Método
                     hessian = T))                  # Calcular a hessiana
    
    if(typeof(op) == 'character')
    { op      <- try(optim(par = par, # Chute inicial
                           fn = log_lindley_geometrica,   # Log-Verossimilhança
                           x = amostra,                   # Amostra
                           control = list(fnscale = -1),
                           method = 'BFGS',        # Método
                           hessian = F))    
    
    }
    
    if(typeof(op) == 'character'){op <- list(par = c(NA, NA), convergence = 99)}
    
    h <- try(solve(op$hessian), T)              # Tenta inverter a hessiana
    if(typeof(h) == 'character') {h <- c(NA, NA, NA, NA)}  # Se não for invertível, ele guarda o erro em character
    # Daí se o tipo for character, h vira um vetor de NA
    
    valores <- c(op$par[1], op$par[2], par[1], par[2], n, -h[1], -h[4], rep(0,8))
    # Valores recebe o que queremos dessa bagaça toda,
    #theta_estimado, rho_estimado, theta_real, rho_real, n, variância_rho, variância_theta, e os zeros serão substituídos fora do for por vício_theta, eqm_theta, erro padrão_theta, probabilidade de cobertura_theta, vício_rho, eqm_rho, erro padrão_rho, probabilidade de cobertura_rhoa
    
    cat('itr:', i, '-' , valores, '\n')
    simulacoes_bfgs[i, ,index_par, index_n] <- valores
    
    }
  }
  
}


simulacoes_bfgs


simulacoes_bfgs[,,9,1]


simulacoes_bfgs[,8,,]  <-  simulacoes_bfgs[,1,,] - simulacoes_bfgs[,3,,]     # Vício de rho
simulacoes_bfgs[,9,,]  <- (simulacoes_bfgs[,1,,] - simulacoes_bfgs[,3,,])^2 # EQM de rho
simulacoes_bfgs[,10,,] <- (simulacoes_bfgs[,6,,])^0.5
simulacoes_bfgs[,11,,] <- ((simulacoes_bfgs[,1,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_bfgs[,6,,])) < simulacoes_bfgs[,3,,]) & ((simulacoes_bfgs[,1,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_bfgs[,6,,])) > simulacoes_bfgs[,3,,]) # Prob cobertura rho

simulacoes_bfgs[,12,,] <-  simulacoes_bfgs[,2,,] - simulacoes_bfgs[,4,,]
simulacoes_bfgs[,13,,]  <- (simulacoes_bfgs[,2,,] - simulacoes_bfgs[,4,,])^2
simulacoes_bfgs[,14,,] <- (simulacoes_bfgs[,7,,])^0.5
simulacoes_bfgs[,15,,] <- ((simulacoes_bfgs[,2,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_bfgs[,7,,])) < simulacoes_bfgs[,4,,]) & ((simulacoes_bfgs[,2,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_bfgs[,7,,])) > simulacoes_bfgs[,4,,]) # Prob cobertura rho

simulacoes_bfgs

diagnostico_bfgs <- simulacoes_bfgs[,c(5:15),,]

save(simulacoes_bfgs,  file = 'simulacoes_bfgs.Rdata')
save(diagnostico_bfgs, file = 'diagnostico_bfgs.Rdata')


# Método CG ---------------------------------------------------------------

library(LambertW)

# Uma tentativa de maximização com optim
optim(par = c(1, 0.3), fn = log_lindley_geometrica, x = rlingley_geom(1000, c(4, 0.1)), control = list(fnscale = -1), method = 'CG', hessian = T)

par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000

simulacoes_cg <- array(c(rep(0,6)), dim=c(N,15,9,10))
simulacoes_cg

dim(simulacoes_cg) # Serão 50 mil linhas, 15 colunas e 90 matrizes
# As dimensões representam, respectivamente, Linha, coluna, dimensão referente a comb dos parâmetros, dimensão do tamanho de amostra

set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{ 
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9)                         # Combinação de parâmetros
  { par <- par_comb[[index_par]]
  amostra <- rlingley_geom(n, par=par)           # Amostra
  op      <- try(optim(par = par,            # Chute inicial
                       fn = log_lindley_geometrica,  # Log-Verossimilhança
                       x = amostra,                  # Amostra
                       control = list(fnscale = -1),
                       method = 'CG',              # Método
                       hessian = T))                  # Calcular a hessiana
  
  if(typeof(op) == 'character')
  { op      <- try(optim(par = par, # Chute inicial
                         fn = log_lindley_geometrica,   # Log-Verossimilhança
                         x = amostra,                   # Amostra
                         control = list(fnscale = -1),
                         method = 'CG',        # Método
                         hessian = F))                  
  }
  
  if(typeof(op) == 'character')
  { valores <- c(NA, NA, par[1], par[2], n, NA, NA,  rep(0,8))
  simulacoes_cg[i, ,index_par, index_n] <- valores
  next}
  
  h <- try(solve(op$hessian))              # Tenta inverter a hessiana
  if(typeof(h) == 'character') {h <- c(NA, NA, NA, NA)}  # Se não for invetível, ele guarda o erro em character
  # Daí se o tipo for character, h vira um vetor de NA
  
  valores <- c(op$par[1], op$par[2], par[1], par[2], n, -h[1], -h[4], rep(0,8))
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
  
  cat('itr:', i, '-' , valores, '\n')
  simulacoes_cg[i, ,index_par, index_n] <- valores
  
  }
  }
  
}

simulacoes_cg


simulacoes_cg[,8,,]  <-  simulacoes_cg[,1,,] - simulacoes_cg[,3,,]     # Vício de rho
simulacoes_cg[,9,,]  <- (simulacoes_cg[,1,,] - simulacoes_cg[,3,,])^2 # EQM de rho
simulacoes_cg[,10,,] <- (simulacoes_cg[,6,,])^0.5
simulacoes_cg[,11,,] <- ((simulacoes_cg[,1,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_cg[,6,,])) < simulacoes_cg[,3,,]) & ((simulacoes_cg[,1,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_cg[,6,,])) > simulacoes_cg[,3,,]) # Prob cobertura rho

simulacoes_cg[,12,,] <-  simulacoes_cg[,2,,] - simulacoes_cg[,4,,]
simulacoes_cg[,13,,]  <- (simulacoes_cg[,2,,] - simulacoes_cg[,4,,])^2
simulacoes_cg[,14,,] <- (simulacoes_cg[,7,,])^0.5
simulacoes_cg[,15,,] <- ((simulacoes_cg[,2,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_cg[,7,,])) < simulacoes_cg[,4,,]) & ((simulacoes_cg[,2,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_cg[,7,,])) > simulacoes_cg[,4,,]) # Prob cobertura rho

simulacoes_cg

diagnostico_cg <- simulacoes_cg[,c(5:15),,]

save(simulacoes_cg,  file = 'simulacoes_cg.Rdata')
save(diagnostico_cg, file = 'diagnostico_cg.Rdata')

# Método L-BFGS-B ---------------------------------------------------------

library(LambertW)

# Uma tentativa de maximização com optim

optim(par = c(1, 0.3), fn = log_lindley_geometrica, x = rlingley_geom(1000, c(4, 0.1)),
      control = list(fnscale = -1), method = 'L-BFGS-B', lower = c(1.00001,0.001), upper = c(Inf,1), hessian = T)

# Me parece meio instável, talvez se reduzir o intervalo
optim(par = c(2, 0.3), fn = log_lindley_geometrica, x = rlingley_geom(1000, c(4, 0.1)),
      control = list(fnscale = -1), method = 'L-BFGS-B', lower = c(0,0.001), upper = c(Inf,1), hessian = T)

par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000

simulacoes_lbfgs <- array(c(rep(0,6)), dim=c(N,15,9,10))
simulacoes_lbfgs


set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{ 
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9)                         # Combinação de parâmetros
  { par <- par_comb[[index_par]]
  amostra <- rlingley_geom(n, par=par)               # Amostra
  op      <- try(optim(par = par,            # Chute inicial
                       fn = log_lindley_geometrica,  # Log-Verossimilhança
                       x = amostra,                  # Amostra
                       control = list(fnscale = -1),
                       method = 'L-BFGS-B',          # Método
                       hessian = T,                  # Calcular a hessiana
                       lower = c(0.000001, 0.000001),
                       upper = c(100, 0.999999)
                       ))                 
  
  if(typeof(op) == 'character')
  { op      <-  try(optim(par = par,            # Chute inicial
                          fn = log_lindley_geometrica,  # Log-Verossimilhança
                          x = amostra,                  # Amostra
                          control = list(fnscale = 1),
                          method = 'L-BFGS-B',          # Método
                          lower = c(0.000001,  0.000001),
                          upper = c(100, 0.999999)
  ))                   
  }
  
  if(typeof(op) == 'character')
  { valores <- c(NA, NA, par[1], par[2], n, NA, NA,  rep(0,8))
    simulacoes_lbfgs[i, ,index_par, index_n] <- valores
    next}
  
  h <- try(solve(op$hessian))              # Tenta inverter a hessiana
  if(typeof(h) == 'character') {h <- c(NA, NA, NA, NA)}  # Se não for invetível, ele guarda o erro em character
  # Daí se o tipo for character, h vira um vetor de NA
  
  valores <- c(op$par[1], op$par[2], par[1], par[2], n, -h[1], -h[4], rep(0,8))
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
  
  cat('itr:', i, '-' , valores, '\n')
  simulacoes_lbfgs[i, ,index_par, index_n] <- valores
  
  }
  }
  
}

simulacoes_lbfgs

simulacoes_lbfgs[,8,,]  <-  simulacoes_lbfgs[,1,,] - simulacoes_lbfgs[,3,,]     # Vício de rho
simulacoes_lbfgs[,9,,]  <- (simulacoes_lbfgs[,1,,] - simulacoes_lbfgs[,3,,])^2 # EQM de rho
simulacoes_lbfgs[,10,,] <- (simulacoes_lbfgs[,6,,])^0.5
simulacoes_lbfgs[,11,,] <- ((simulacoes_lbfgs[,1,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_lbfgs[,6,,])) < simulacoes_lbfgs[,3,,]) & ((simulacoes_lbfgs[,1,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_lbfgs[,6,,])) > simulacoes_lbfgs[,3,,]) # Prob cobertura rho

simulacoes_lbfgs[,12,,] <-  simulacoes_lbfgs[,2,,] - simulacoes_lbfgs[,4,,]
simulacoes_lbfgs[,13,,]  <- (simulacoes_lbfgs[,2,,] - simulacoes_lbfgs[,4,,])^2
simulacoes_lbfgs[,14,,] <- (simulacoes_lbfgs[,7,,])^0.5
simulacoes_lbfgs[,15,,] <- ((simulacoes_lbfgs[,2,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_lbfgs[,7,,])) < simulacoes_lbfgs[,4,,]) & ((simulacoes_lbfgs[,2,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_lbfgs[,7,,])) > simulacoes_lbfgs[,4,,]) # Prob cobertura rho

simulacoes_lbfgs

diagnostico_lbfgs <- simulacoes_lbfgs[,c(5:15),,]

save(simulacoes_lbfgs,  file = 'simulacoes_lbfgs.Rdata')
save(diagnostico_lbfgs, file = 'diagnostico_lbfgs.Rdata')


# SANN --------------------------------------------------------------------

library(LambertW)
amostra <- rlingley_geom(10000, c(5, 0.8))

optim(par = c(1, 0.3), fn = log_lindley_geometrica, x = amostra, control = list(fnscale = -1))
#Restrições: par[1] > 0 e 0 < par[2] < 1


par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8),   # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),   # Daí são as combinações
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000                                                # O N de simulações pedido é 50000

simulacoes_sann <- array(c(rep(0,6)), dim=c(N,15,9,10))  # Esse array vai guardar os resultados
simulacoes_sann                                         # Basicamente são 90 matrizes dentro de um array

dim(simulacoes_sann) # Serão 50 mil linhas, 15 colunas e 90 matrizes

# As dimensões representam, respectivamente, Linha, coluna, dimensão referente a comb dos parâmetros, dimensão do tamanho de amostra
set.seed(9999)
for (i in 1:N) # Número de simulações
{ 
  for (index_n in 1:10) # Tamanho da amostra
  {n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9) # Combinação de parâmetros
  { par <- par_comb[[index_par]]
  amostra <- rlingley_geom(n, par=par)     # Amostra
  op      <- try(optim(par = par, # Chute inicial
                       fn = log_lindley_geometrica,   # Log-Verossimilhança
                       x = amostra,                   # Amostra
                       control = list(fnscale = -1),
                       method = 'SANN',        # Método
                       hessian = T), T)                  # Calcular a hessiana
  
  if(typeof(op) == 'character')
  { op      <- try(optim(par = par, # Chute inicial
                         fn = log_lindley_geometrica,   # Log-Verossimilhança
                         x = amostra,                   # Amostra
                         control = list(fnscale = -1),
                         method = 'SANN',        # Método
                         hessian = F)) 
  
  }
  
  if(typeof(op) == 'character')
  { valores <- c(NA, NA, par[1], par[2], n,  NA, NA, rep(0,8))
  simulacoes_cg[i, ,index_par, index_n] <- valores
  next}
  
  
  
  h <- try(solve(op$hessian))              # Tenta inverter a hessiana
  if(typeof(h) == 'character') {h <- c(NA, NA, NA, NA)}  # Se não for invetível, ele guarda o erro em character
  # Daí se o tipo for character, h vira um vetor de NA
  
  valores <- c(op$par[1], op$par[2], par[1], par[2], n, -h[1], -h[4], rep(0,8))
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, variância_rho, variância_theta, e os zeros serão substituídos fora do for por vício_theta, eqm_theta, erro padrão_theta, probabilidade de cobertura_theta, vício_rho, eqm_rho, erro padrão_rho, probabilidade de cobertura_rho
  
  cat('itr:', i, '-' , valores, '\n')    # Inútil, é só pra vc ficar vendo oq ta acontecendo
  
  simulacoes_sann[i, ,index_par, index_n] <- valores  # Guarda na tabela
  
  }
  }
  
}

simulacoes_sann #theta_estimado, rho_estimado, theta_real, rho_real, n, variância_rho, variância_theta, e os zeros serão substituídos fora do for por vício_theta, eqm_theta, erro padrão_theta, probabilidade de cobertura_theta, vício_rho, eqm_rho, erro padrão_rho, probabilidade de cobertura_rho


# As dimensões representam, respectivamente, Linha, coluna, dimensão referente a comb dos parâmetros, dimensão do tamanho de amostra
simulacoes_sann[,,9,1]


simulacoes_sann[,8,,]  <-  simulacoes_sann[,1,,] - simulacoes_sann[,3,,]     # Vício de rho
simulacoes_sann[,9,,]  <- (simulacoes_sann[,1,,] - simulacoes_sann[,3,,])^2 # EQM de rho
simulacoes_sann[,10,,] <- (simulacoes_sann[,6,,])^0.5
simulacoes_sann[,11,,] <- ((simulacoes_sann[,1,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_sann[,6,,])) < simulacoes_sann[,3,,]) & ((simulacoes_sann[,1,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_sann[,6,,]) > simulacoes_sann[,3,,])) # Prob cobertura rho

simulacoes_sann[,12,,] <-  simulacoes_sann[,2,,] - simulacoes_sann[,4,,]
simulacoes_sann[,13,,]  <- (simulacoes_sann[,2,,] - simulacoes_sann[,4,,])^2
simulacoes_sann[,14,,] <- (simulacoes_sann[,7,,])^0.5
simulacoes_sann[,15,,] <- ((simulacoes_sann[,1,,] - qnorm(1 - 0.05/2) * sqrt(simulacoes_sann[,7,,])) < simulacoes_sann[,4,,]) & ((simulacoes_sann[,1,,] + qnorm(1 - 0.05/2) * sqrt(simulacoes_sann[,7,,])) > simulacoes_sann[,4,,]) # Prob cobertura rho

simulacoes_sann


# Os Resultados -----------------------------------------------------------



# Resultados nelder
for(i in 1:10){
  vicio_theta <- mean(simulacoes_nelder[,8,,i],na.rm=T)
  eqm_theta <- mean(simulacoes_nelder[,9,,i],na.rm=T)
  erro_theta <- mean(simulacoes_nelder[,10,,i],na.rm=T)
  prob_theta <- mean(simulacoes_nelder[,11,,i],na.rm=T)
  var_theta  <- mean(simulacoes_nelder[,6,,i], na.rm = T)
  tam_theta  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_nelder[,6,,i]), na.rm = T)
  
  vicio_rho <- mean(simulacoes_nelder[,12,,i],na.rm=T)
  eqm_rho <- mean(simulacoes_nelder[,13,,i],na.rm=T)
  erro_rho <- mean(simulacoes_nelder[,14,,i],na.rm=T)
  prob_rho <- mean(simulacoes_nelder[,15,,i],na.rm=T)
  var_rho  <- mean(simulacoes_nelder[,7,,i], na.rm = T)
  tam_rho  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_nelder[,7,,i]), na.rm = T)
  
  valores <- c(var_theta, vicio_theta, eqm_theta, erro_theta, prob_theta, tam_theta, var_rho, vicio_rho, eqm_rho, erro_rho, prob_rho, tam_rho, i*10, 'Nelder')
  if(i == 1) 
  {
    resultados_nelder <- data.frame('var_t' = var_theta,'vicio_t'= vicio_theta, 'eqm_t'= eqm_theta,
                                    'erro_t' = erro_theta,  'prob_t' = prob_theta, 'tam_t'= tam_theta,
                                    'var_r' = var_rho,'vicio_r' = vicio_rho, 'eqm_r' = eqm_rho,
                                    'erro_r' = erro_rho, 'prob_r' = prob_rho,'tam_r'= tam_rho,
                                    'n' = i*10, 'simulacao' = 'Nelder')
  }
  else{
    resultados_nelder <-rbind(resultados_nelder, valores)
  }
  
  
}
resultados_nelder

## RESULTADOS BFGS

for(i in 1:10){
  vicio_theta <- mean(simulacoes_bfgs[,8,,i],na.rm=T)
  eqm_theta <- mean(simulacoes_bfgs[,9,,i],na.rm=T)
  erro_theta <- mean(simulacoes_bfgs[,10,,i],na.rm=T)
  prob_theta <- mean(simulacoes_bfgs[,11,,i],na.rm=T)
  var_theta  <- mean(simulacoes_bfgs[,6,,i], na.rm = T)
  tam_theta  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_bfgs[,6,,i]), na.rm = T)
  
  vicio_rho <- mean(simulacoes_bfgs[,12,,i],na.rm=T)
  eqm_rho <- mean(simulacoes_bfgs[,13,,i],na.rm=T)
  erro_rho <- mean(simulacoes_bfgs[,14,,i],na.rm=T)
  prob_rho <- mean(simulacoes_bfgs[,15,,i],na.rm=T)
  var_rho  <- mean(simulacoes_bfgs[,7,,i], na.rm = T)
  tam_rho  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_bfgs[,7,,i]), na.rm = T)
  
  valores <- c(var_theta, vicio_theta, eqm_theta, erro_theta, prob_theta, tam_theta, var_rho, vicio_rho, eqm_rho, erro_rho, prob_rho, tam_rho, i*10, 'BFGS')
  if(i == 1) 
  {
    resultados_bfgs <- data.frame('var_t' = var_theta,'vicio_t'= vicio_theta, 'eqm_t'= eqm_theta,
                                  'erro_t' = erro_theta,  'prob_t' = prob_theta, 'tam_t'= tam_theta,
                                  'var_r' = var_rho,'vicio_r' = vicio_rho, 'eqm_r' = eqm_rho,
                                  'erro_r' = erro_rho, 'prob_r' = prob_rho,'tam_r'= tam_rho,
                                  'n' = i*10, 'simulacao' = 'BFGS')
  }
  else{
    resultados_bfgs <-rbind(resultados_bfgs, valores)
  }
  
  
}
resultados_bfgs

## RESULTADOS CG

for(i in 1:10){
  vicio_theta <- mean(simulacoes_cg[,8,,i],na.rm=T)
  eqm_theta <- mean(simulacoes_cg[,9,,i],na.rm=T)
  erro_theta <- mean(simulacoes_cg[,10,,i],na.rm=T)
  prob_theta <- mean(simulacoes_cg[,11,,i],na.rm=T)
  var_theta  <- mean(simulacoes_cg[,6,,i], na.rm = T)
  tam_theta  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_cg[,6,,i]), na.rm = T)
  
  vicio_rho <- mean(simulacoes_cg[,12,,i],na.rm=T)
  eqm_rho <- mean(simulacoes_cg[,13,,i],na.rm=T)
  erro_rho <- mean(simulacoes_cg[,14,,i],na.rm=T)
  prob_rho <- mean(simulacoes_cg[,15,,i],na.rm=T)
  var_rho  <- mean(simulacoes_cg[,7,,i], na.rm = T)
  tam_rho  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_cg[,7,,i]), na.rm = T)
  
  valores <- c(var_theta, vicio_theta, eqm_theta, erro_theta, prob_theta, tam_theta, var_rho, vicio_rho, eqm_rho, erro_rho, prob_rho, tam_rho, i*10, 'CG')
  if(i == 1) 
  {
    resultados_cg <- data.frame('var_t' = var_theta,'vicio_t'= vicio_theta, 'eqm_t'= eqm_theta,
                                'erro_t' = erro_theta,  'prob_t' = prob_theta, 'tam_t'= tam_theta,
                                'var_r' = var_rho,'vicio_r' = vicio_rho, 'eqm_r' = eqm_rho,
                                'erro_r' = erro_rho, 'prob_r' = prob_rho,'tam_r'= tam_rho,
                                'n' = i*10, 'simulacao' = 'CG')
  }
  else{
    resultados_cg <-rbind(resultados_cg, valores)
  }
  
  
}
resultados_cg

## RESULTADOS LBFGS

for(i in 1:10){
  vicio_theta <- mean(simulacoes_lbfgs[,8,,i],na.rm=T)
  eqm_theta <- mean(simulacoes_lbfgs[,9,,i],na.rm=T)
  erro_theta <- mean(simulacoes_lbfgs[,10,,i],na.rm=T)
  prob_theta <- mean(simulacoes_lbfgs[,11,,i],na.rm=T)
  var_theta  <- mean(simulacoes_lbfgs[,6,,i], na.rm = T)
  tam_theta  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_lbfgs[,6,,i]), na.rm = T)
  
  vicio_rho <- mean(simulacoes_lbfgs[,12,,i],na.rm=T)
  eqm_rho <- mean(simulacoes_lbfgs[,13,,i],na.rm=T)
  erro_rho <- mean(simulacoes_lbfgs[,14,,i],na.rm=T)
  prob_rho <- mean(simulacoes_lbfgs[,15,,i],na.rm=T)
  var_rho  <- mean(simulacoes_lbfgs[,7,,i], na.rm = T)
  tam_rho  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_lbfgs[,7,,i]), na.rm = T)
  
  valores <- c(var_theta, vicio_theta, eqm_theta, erro_theta, prob_theta, tam_theta, var_rho, vicio_rho, eqm_rho, erro_rho, prob_rho, tam_rho, i*10, 'L-BFGS-B')
  if(i == 1) 
  {
    resultados_lbfgs <- data.frame('var_t' = var_theta,'vicio_t'= vicio_theta, 'eqm_t'= eqm_theta,
                                   'erro_t' = erro_theta,  'prob_t' = prob_theta, 'tam_t'= tam_theta,
                                   'var_r' = var_rho,'vicio_r' = vicio_rho, 'eqm_r' = eqm_rho,
                                   'erro_r' = erro_rho, 'prob_r' = prob_rho,'tam_r'= tam_rho,
                                   'n' = i*10, 'simulacao' = 'L-BFGS-B')
  }
  else{
    resultados_lbfgs <-rbind(resultados_lbfgs, valores)
  }
  
  
}
resultados_lbfgs

## RESULTADOS SANN

for(i in 1:10){
  vicio_theta <- mean(simulacoes_sann[,8,,i],na.rm=T)
  eqm_theta <- mean(simulacoes_sann[,9,,i],na.rm=T)
  erro_theta <- mean(simulacoes_sann[,10,,i],na.rm=T)
  prob_theta <- mean(simulacoes_sann[,11,,i],na.rm=T)
  var_theta  <- mean(simulacoes_sann[,6,,i], na.rm = T)
  tam_theta  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_sann[,6,,i]), na.rm = T)
  
  vicio_rho <- mean(simulacoes_sann[,12,,i],na.rm=T)
  eqm_rho <- mean(simulacoes_sann[,13,,i],na.rm=T)
  erro_rho <- mean(simulacoes_sann[,14,,i],na.rm=T)
  prob_rho <- mean(simulacoes_sann[,15,,i],na.rm=T)
  var_rho  <- mean(simulacoes_sann[,7,,i], na.rm = T)
  tam_rho  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_sann[,7,,i]), na.rm=T)
  
  valores <- c(var_theta, vicio_theta, eqm_theta, erro_theta, prob_theta, tam_theta, var_rho, vicio_rho, eqm_rho, erro_rho, prob_rho, tam_rho, i*10, 'SANN')
  if(i == 1) 
  {
    resultados_sann <- data.frame('var_t' = var_theta,'vicio_t'= vicio_theta, 'eqm_t'= eqm_theta,
                                  'erro_t' = erro_theta,  'prob_t' = prob_theta, 'tam_t'= tam_theta,
                                  'var_r' = var_rho,'vicio_r' = vicio_rho, 'eqm_r' = eqm_rho,
                                  'erro_r' = erro_rho, 'prob_r' = prob_rho,'tam_r'= tam_rho,
                                  'n' = i*10, 'simulacao' = 'SANN')
  }
  else{
    resultados_sann <-rbind(resultados_sann, valores)
  }
  
  
}
resultados_sann

## RESULTADOS NEWTON-RAPHSON

for(i in 1:10){
  vicio_theta <- mean(simulacoes_newton[,8,,i],na.rm=T)
  eqm_theta <- mean(simulacoes_newton[,9,,i],na.rm=T)
  erro_theta <- mean(simulacoes_newton[,10,,i],na.rm=T)
  prob_theta <- mean(simulacoes_newton[,11,,i],na.rm=T)
  var_theta  <- mean(simulacoes_newton[,6,,i], na.rm = T)
  tam_theta  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_newton[,6,,i]), na.rm = T)
  
  vicio_rho <- mean(simulacoes_newton[,12,,i],na.rm=T)
  eqm_rho <- mean(simulacoes_newton[,13,,i],na.rm=T)
  erro_rho <- mean(simulacoes_newton[,14,,i],na.rm=T)
  prob_rho <- mean(simulacoes_newton[,15,,i],na.rm=T)
  var_rho  <- mean(simulacoes_newton[,7,,i], na.rm = T)
  tam_rho  <- 2*qnorm(1 - 0.05/2)*mean(sqrt(simulacoes_newton[,7,,i]), ra.rm = T)
  
  valores <- c(var_theta, vicio_theta, eqm_theta, erro_theta, prob_theta, tam_theta, var_rho, vicio_rho, eqm_rho, erro_rho, prob_rho, tam_rho, i*10, 'NR')
  if(i == 1) 
  {
    resultados_newton <- data.frame('var_t' = var_theta,'vicio_t'= vicio_theta, 'eqm_t'= eqm_theta,
                                    'erro_t' = erro_theta,  'prob_t' = prob_theta, 'tam_t'= tam_theta,
                                    'var_r' = var_rho,'vicio_r' = vicio_rho, 'eqm_r' = eqm_rho,
                                    'erro_r' = erro_rho, 'prob_r' = prob_rho,'tam_r'= tam_rho, 
                                    'n' = i*10, 'simulacao' = 'NR')
  }
  else{
    resultados_newton <-rbind(resultados_newton, valores)
  }
  
  
}
resultados_newton



# Análise dos resultados --------------------------------------------------

resultados_nelder
resultados_bfgs
resultados_cg
resultados_lbfgs
resultados_sann
resultados_newton



# UNIVARIADA

library(ggplot2)

uni_lineplot <- function(var, banco, simulacao,  ylab = '', xlab = 'Tamanho da amostra')
{
  if(var == 'prob_t' || var == 'prob_r')
  {
    grafico <- ggplot(banco,aes(x= n, y= .data[[var]])) +
      geom_line()+ 
      geom_point()+ 
      geom_hline(yintercept=0.95, col = 'red', linetype='dashed')+
      theme_minimal()+
      xlab(xlab)+
      ylab(ylab)
  }
  else{
    grafico <- ggplot(banco,aes(x= n, y= .data[[var]])) +
      geom_line()+ 
      geom_point()+ 
      geom_hline(yintercept=0, col = 'red', linetype='dashed')+
      theme_minimal()+
      xlab(xlab)+
      ylab(ylab)}
  
  print(grafico)
  arquivo <- paste0('img/', simulacao, '_', var, '.jpg')
  ggsave(arquivo)
}

uni_lineplot('vicio_t', resultados_cg, 'cg')
uni_lineplot('eqm_t', resultados_cg, 'cg')
uni_lineplot('prob_t', resultados_cg, 'cg')

# COMPARATIVO -------------------------------

resultados <- do.call(rbind,list(resultados_nelder,
                                 resultados_bfgs,
                                 resultados_cg,
                                 resultados_lbfgs,
                                 resultados_newton
))

resultados$n <- factor(resultados$n, levels=c(seq(10, 100,10)))
save(resultados, file = 'resultados.RData')

library(ggplot2)
library(gridExtra)
formato <- theme(                                                       
  plot.title = element_text(size = 14, hjust = 0.5),
  axis.title.y = element_text(size = 12, vjust = 0.5, angle= 0),
  axis.title.x = element_text(size = 12, vjust = -0.2),
  axis.text.y = element_text(size = 10),
  axis.text.x = element_text(size = 10))

paleta <- c("#662F00", "#996035", "#CC9B7A", '#D8AF97')

par(mfrow = c(1,1))

# Vício para theta
ggplot(resultados, aes(x = n, y = as.numeric(vicio_t), colour = simulacao, group = simulacao)) +
  geom_line(linewidth =1.3, alpha = 0.7)+
  geom_point(size=2)+
  xlab('Tamanho da amostra')+
  ylab('Vício')+
  ggtitle('Gráfico do vício para theta por simulação')+
  theme_minimal()+
  formato +
  guides(colour = guide_legend(title = "Simulação"))+
  scale_colour_manual(values=paleta)+
  geom_hline(yintercept=0, col = 'red', linetype = 'dashed')

ggsave('img/Vicio_theta.jpg')

# Vício para rho
ggplot(resultados, aes(x = n, y = as.numeric(vicio_r), colour = simulacao, group = simulacao)) +
  geom_line(linewidth =1.3, alpha = 0.7)+
  geom_point(size=2)+
  xlab('Tamanho da amostra')+
  ylab('Vício')+
  ggtitle('Gráfico do vício para rho por simulação')+
  theme_minimal()+
  formato +
  guides(colour = guide_legend(title = "Simulação"))+
  scale_colour_manual(values=paleta)+
  geom_hline(yintercept=0, col = 'red', linetype = 'dashed')


ggsave('img/Vicio_rho.jpg')

# eqm para theta
ggplot(resultados, aes(x = n, y = as.numeric(eqm_t), colour = simulacao, group = simulacao)) +
  geom_line(linewidth =1.3, alpha = 0.7)+
  geom_point(size=2)+
  xlab('Tamanho da amostra')+
  ylab('EQM')+
  ggtitle('Gráfico do EQM para theta por simulação')+
  theme_minimal()+
  formato +
  guides(colour = guide_legend(title = "Simulação"))+
  scale_colour_manual(values=paleta)+
  geom_hline(yintercept=0, col = 'red', linetype = 'dashed')

ggsave('img/eqm_theta.jpg')

# eqm para rho
ggplot(resultados, aes(x = n, y = as.numeric(eqm_r), colour = simulacao, group = simulacao)) +
  geom_line(linewidth =1.3, alpha = 0.7)+
  geom_point(size=2)+
  xlab('Tamanho da amostra')+
  ylab('EQM')+
  ggtitle('Gráfico do EQM para rho por simulação')+
  theme_minimal()+
  formato +
  guides(colour = guide_legend(title = "Simulação"))+
  scale_colour_manual(values=paleta)+
  geom_hline(yintercept=0, col = 'red', linetype = 'dashed')

ggsave('img/eqm_rho.jpg')

# erro para theta
ggplot(resultados, aes(x = n, y = as.numeric(erro_t), colour = simulacao, group = simulacao)) +
  geom_line(linewidth =1, alpha = 0.7)+
  geom_point(size=2)+
  xlab('Tamanho da amostra')+
  ylab('Erro')+
  ggtitle('Gráfico do erro padrão para theta por simulação')+
  theme_minimal()+
  formato +
  guides(colour = guide_legend(title = "Simulação"))+
  scale_colour_manual(values=paleta)+
  geom_hline(yintercept=0, col = 'red', linetype = 'dashed')

ggsave('img/erro_theta.jpg')

# erro para rho
ggplot(resultados, aes(x = n, y = as.numeric(erro_r), colour = simulacao, group = simulacao)) +
  geom_line(linewidth =1, alpha = 0.7)+
  geom_point(size=2)+
  xlab('Tamanho da amostra')+
  ylab('Erro')+
  ggtitle('Gráfico do erro padrão para rho por simulação')+
  theme_minimal()+
  formato +
  guides(colour = guide_legend(title = "Simulação"))+
  scale_colour_manual(values=paleta)+
  geom_hline(yintercept=0, col = 'red', linetype = 'dashed')


ggsave('img/erro_rho.jpg') 

# probabilidade de cobertura para theta
ggplot(resultados, aes(x = n, y = as.numeric(prob_t), colour = simulacao, group = simulacao)) +
  geom_line(linewidth =1, alpha = 0.7)+
  geom_point(size=2)+
  xlab('Tamanho da amostra')+
  ylab('Probabilidade')+
  ggtitle('Gráfico da probabilidade de cobertura para theta por simulação')+
  theme_minimal()+
  formato +
  guides(colour = guide_legend(title = "Simulação"))+
  scale_colour_manual(values=paleta)+
  geom_hline(yintercept=0.95, col = 'red', linetype = 'dashed')

ggsave('img/prob_theta.jpg')

# probabilidade de cobertura para rho
ggplot(resultados, aes(x = n, y = as.numeric(prob_r), colour = simulacao, group = simulacao)) +
  geom_line(linewidth =1, alpha = 0.7)+
  geom_point(size=2)+
  xlab('Tamanho da amostra')+
  ylab('Probabilidade')+
  ggtitle('Gráfico da probabilidade de cobertura para rho por simulação')+
  theme_minimal()+
  formato +
  guides(colour = guide_legend(title = "Simulação"))+
  scale_colour_manual(values=paleta)+
  geom_hline(yintercept=0.95, col = 'red', linetype = 'dashed')

ggsave('img/prob_rho.jpg')

# tamanho de cobertura para theta
ggplot(resultados, aes(x = n, y = as.numeric(tam_t), colour = simulacao, group = simulacao)) +
  geom_line(linewidth =1, alpha = 0.7)+
  geom_point(size=2)+
  xlab('Tamanho da amostra')+
  ylab('Cobertura')+
  ggtitle('Gráfico do tamanho da cobertura para theta por simulação')+
  theme_minimal()+
  formato +
  guides(colour = guide_legend(title = "Simulação"))+
  scale_colour_manual(values=paleta)+
  geom_hline(yintercept=0, col = 'red', linetype = 'dashed')

ggsave('img/tamanho_theta.jpg')

# tamanho de cobertura para rho
ggplot(resultados, aes(x = n, y = as.numeric(tam_r), colour = simulacao, group = simulacao)) +
  geom_line(linewidth =1, alpha = 0.7)+
  geom_point(size=2)+
  xlab('Tamanho da amostra')+
  ylab('Cobertura')+
  ggtitle('Gráfico do tamanho da cobertura para rho por simulação')+
  theme_minimal()+
  formato +
  guides(colour = guide_legend(title = "Simulação"))+
  scale_colour_manual(values=paleta)+
  geom_hline(yintercept=0, col = 'red', linetype = 'dashed')

ggsave('img/tamanho_rho.jpg')


