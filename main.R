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

rlingley_geom <- function(n, par)
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
  
  n       <- length(x)
  
  dtheta1 <- (((2*x^2+2*x)*(par[1]^2)+(2*x^2+4*x)*par[1])*par[2])
  dtheta2 <- (((x+1)*par[1]^2+(x+2)*par[1]+1)*par[2]+(-(par[1]^2)-2*par[1]-1)*exp(x*par[1]))
  dtheta  <- 2*n/par[1] - n/(par[1] + 1) - sum(x) + sum(dtheta1/dtheta2)
  
  drho1   <- ((x*par[1])/(par[1]+1)+1)
  drho2   <- exp(-(x*par[1]))
  drho    <- -n/(1- par[2]) + sum((2*drho1*drho2)/(1 - drho1*drho2*par[2]))
  
  vetor   <- matrix(nrow = 2, c(dtheta, drho))
 
  return(vetor)
}

set.seed(123)
U(rlingley_geom(1000, 3, 0.5), c(3, 0.5))


### Implementação da matriz de Informação de Fisher - Jacobiano

J <- function(x, par)
{
  
  jacobiano      <- matrix(nrow = 2, ncol = 2)
  jacobiano[1,1] <- dtheta2(x, par)
  jacobiano[1,2] <- dthetarho(x, par)
  jacobiano[2,1] <- jacobiano[1,2]
  jacobiano[2,2] <- drho2(x, par)
  
  return(-jacobiano)
  
}

J(rlingley_geom(1000, 0.3, 0.5), c(0.3, 0.5))


### Função Escore de Fisher 

escore_f <- function(x, par, itr, erro = 0.00001)     # Comecei a implementar uma função, mas tá explodindo
{                                                     # daí to testando fora da função, logo a baixo
  par <- matrix(par, nrow = 2)
  
  i <- 1
  
  while(i < itr)
    
  {
    i <- i + 1
    
    par_novo <- par - solve(J(x, par)) %*% U(x, par)
    
    if(max(abs(U(x, par) - U(x, par))) < erro) break # Critério de parada
    
    par <- par_novo
    
    a <- cat('Iteração:', i, 'shape:', par_novo[1], 'scale:', par_novo[2], "\n")
    
  }
  return(a)
  
}



par <- matrix(c(6, 0.3), nrow = 2)                # Os valores que o Newton-raphson precisa
x <-  rlingley_geom(1000, c(4, 0.1))
itr <- 5000
erro <- 0.00001
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

simulacoes_nelder <- array(c(rep(0,6)), dim=c(N,8,9,10))  # Esse array vai guardar os resultados
simulacoes_nelder                                         # Basicamente são 90 matrizes dentro de um array
 
dim(simulacoes_nelder) # Serão 50 mil linhas, 8 colunas e 90 matrizes

# As dimensões representam, respectivamente, Linha, coluna, dimensão referente a comb dos parâmetros, dimensão do tamanho de amostra


<<<<<<< HEAD

=======
set.seed(9999)
index_par <- 1


#BREAKPOINT
>>>>>>> eb4dc0105d987a8ae7aebdd973062335bae5eebb
for (i in 1:N) # Número de simulações
{ set.seed(9999)
  for (index_n in 1:10) # Tamanho da amostra
  {n <- seq(10, 100, 10)[index_n]
  
    for (index_par in 1:9) # Combinação de parâmetros
    { par <- par_comb[[index_par]]
      amostra <- rlingley_geom(n, par=par)     # Amostra
      op      <- try(optim(par = c(1.5, 0.3), # Chute inicial
                fn = log_lindley_geometrica,   # Log-Verossimilhança
                x = amostra,                   # Amostra
                control = list(fnscale = -1),
                method = 'Nelder-Mead',        # Método
                hessian = T), T)                  # Calcular a hessiana
      
      if(typeof(op) == 'character')
      { op      <- try(optim(par = c(1.5, 0.3), # Chute inicial
                             fn = log_lindley_geometrica,   # Log-Verossimilhança
                             x = amostra,                   # Amostra
                             control = list(fnscale = -1),
                             method = 'Nelder-Mead',        # Método
                             hessian = F)) 
      
      }
      
      h <- try(solve(op$hessian))              # Tenta inverter a hessiana
      if(typeof(h) == 'character') {h <- c(NA, NA, NA, NA)}  # Se não for invetível, ele guarda o erro em character
                                                             # Daí se o tipo for character, h vira um vetor de NA
      
      valores <- c(op$par[1], op$par[2], par[1], par[2], n, op$convergence, h[1], h[4])
      # Valores recebe o que queremos dessa bagaça toda,
      # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
      
      cat('itr:', i, '-' , valores, '\n')    # Inútil, é só pra vc ficar vendo oq ta acontecendo
      
      simulacoes_nelder[i, ,index_par, index_n] <- valores  # Guarda na tabela
      
    }
  }
  
}

simulacoes_nelder # é um array que recebe por coluna - theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta

dim(simulacoes_nelder)

# As dimensões representam, respectivamente, Linha, coluna, dimensão referente a comb dos parâmetros, dimensão do tamanho de amostra
simulacoes_nelder[,,9,1]


simulacoes_nelder[,1,,] - simulacoes_nelder[,3,,]  # dá pra fazer assim vicio
(simulacoes_nelder[,1,,] - simulacoes_nelder[,3,,])^2 # EQM

((simulacoes_nelder[,1,,] - qnorm(1 - 0.05/2) * simulacoes_nelder[,7,,]) < simulacoes_nelder[,3,,]) 
& ((simulacoes_nelder[,1,,] + qnorm(1 - 0.05/2) * simulacoes_nelder[,7,,]) > simulacoes_nelder[,3,,]) # IC







# BFGS --------------------------------------------------------------------

library(LambertW)

# Uma tentativa de maximização com optim
optim(par = c(1, 0.3), fn = log_lindley_geometrica, x = rlingley_geom(1000, c(4, 0.1)), control = list(fnscale = -1), method = 'BFGS', hessian = T)
 
par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000

simulacoes_bfgs <- array(c(rep(0,6)), dim=c(N,8,9,10))
simulacoes_bfgs

dim(simulacoes_bfgs) # Serão 50 mil linhas, 8 colunas e 90 matrizes
# As dimensões representam, respectivamente, Linha, coluna, dimensão referente a comb dos parâmetros, dimensão do tamanho de amostra



for (i in 1:N)                                  # Número de simulações
{ set.seed(9999)
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
    for (index_par in 1:9)                         # Combinação de parâmetros
    { par <- par_comb[[index_par]]
    amostra <- rlingley_geom(n, par=par)           # Amostra
    op      <- try(optim(par = c(1.5, 0.3),            # Chute inicial
                     fn = log_lindley_geometrica,  # Log-Verossimilhança
                     x = amostra,                  # Amostra
                     control = list(fnscale = -1),
                     method = 'BFGS',              # Método
                     hessian = T))                  # Calcular a hessiana
    
    if(typeof(op) == 'character')
    { op      <- try(optim(par = c(1.5, 0.3), # Chute inicial
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
    
    valores <- c(op$par[1], op$par[2], par[1], par[2], n, op$convergence, h[1], h[4])
    # Valores recebe o que queremos dessa bagaça toda,
    # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
    
    cat('itr:', i, '-' , valores, '\n')
    simulacoes_bfgs[i, ,index_par, index_n] <- valores
    
    }
  }
  
}
a <- try(solve(asdas))
a
a <- try(log('a'), T) # Só pra entender como try funciona
a

simulacoes_bfgs


# Método CG ---------------------------------------------------------------

library(LambertW)

# Uma tentativa de maximização com optim
optim(par = c(1, 0.3), fn = log_lindley_geometrica, x = rlingley_geom(1000, c(4, 0.1)), control = list(fnscale = -1), method = 'CG', hessian = T)

par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000

simulacoes_cg <- array(c(rep(0,6)), dim=c(N,8,9,10))
simulacoes_cg

dim(simulacoes_cg) # Serão 50 mil linhas, 8 colunas e 90 matrizes
# As dimensões representam, respectivamente, Linha, coluna, dimensão referente a comb dos parâmetros, dimensão do tamanho de amostra

set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{ 
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9)                         # Combinação de parâmetros
  { par <- par_comb[[index_par]]
  amostra <- rlingley_geom(n, par=par)           # Amostra
  op      <- try(optim(par = c(1.5, 0.3),            # Chute inicial
                       fn = log_lindley_geometrica,  # Log-Verossimilhança
                       x = amostra,                  # Amostra
                       control = list(fnscale = -1),
                       method = 'CG',              # Método
                       hessian = T))                  # Calcular a hessiana
  
  if(typeof(op) == 'character')
  { op      <- try(optim(par = c(1.5, 0.3), # Chute inicial
                         fn = log_lindley_geometrica,   # Log-Verossimilhança
                         x = amostra,                   # Amostra
                         control = list(fnscale = -1),
                         method = 'CG',        # Método
                         hessian = F))                  
  }
  
  if(typeof(op) == 'character')
  { valores <- c(NA, NA, par[1], par[2], n, 99, NA, NA)
  simulacoes_cg[i, ,index_par, index_n] <- valores
  next}
  
  h <- try(solve(op$hessian))              # Tenta inverter a hessiana
  if(typeof(h) == 'character') {h <- c(NA, NA, NA, NA)}  # Se não for invetível, ele guarda o erro em character
  # Daí se o tipo for character, h vira um vetor de NA
  
  valores <- c(op$par[1], op$par[2], par[1], par[2], n, op$convergence, h[1], h[4])
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
  
  cat('itr:', i, '-' , valores, '\n')
  simulacoes_cg[i, ,index_par, index_n] <- valores
  
  }
  }
  
}

simulacoes_cg


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

simulacoes_lbfgs <- array(c(rep(0,6)), dim=c(N,8,9,10))
simulacoes_lbfgs


set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{ 
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9)                         # Combinação de parâmetros
  { par <- par_comb[[index_par]]
  amostra <- rlingley_geom(n, par=par)               # Amostra
  op      <- try(optim(par = c(1.5, 0.3),            # Chute inicial
                       fn = log_lindley_geometrica,  # Log-Verossimilhança
                       x = amostra,                  # Amostra
                       control = list(fnscale = -1),
                       method = 'L-BFGS-B',          # Método
                       hessian = T,                  # Calcular a hessiana
                       lower = c(0.000001, 0.000001),
                       upper = c(100, 0.999999)
                       ))                 
  
  if(typeof(op) == 'character')
  { op      <-  try(optim(par = c(1.5, 0.3),            # Chute inicial
                          fn = log_lindley_geometrica,  # Log-Verossimilhança
                          x = amostra,                  # Amostra
                          control = list(fnscale = 1),
                          method = 'L-BFGS-B',          # Método
                          lower = c(0.000001,  0.000001),
                          upper = c(100, 0.999999)
  ))                   
  }
  
  if(typeof(op) == 'character')
  { valores <- c(NA, NA, par[1], par[2], n, 99, NA, NA)
    simulacoes_lbfgs[i, ,index_par, index_n] <- valores
    next}
  
  h <- try(solve(op$hessian))              # Tenta inverter a hessiana
  if(typeof(h) == 'character') {h <- c(NA, NA, NA, NA)}  # Se não for invetível, ele guarda o erro em character
  # Daí se o tipo for character, h vira um vetor de NA
  
  valores <- c(op$par[1], op$par[2], par[1], par[2], n, op$convergence,  h[1], h[4])
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
  
  cat('itr:', i, '-' , valores, '\n')
  simulacoes_lbfgs[i, ,index_par, index_n] <- valores
  
  }
  }
  
}

simulacoes_lbfgs



# Método SANN -------------------------------------------------------------

#install.packages("maxLik")
library(maxLik)

# Esse é o pacote que o professor falou que era mais otimizado pro SANN

# https://cran.r-project.org/web/packages/maxLik/maxLik.pdf                 pag 17

help("maxSANN")

# Pra função do maxlik precisa alterar a LL
log_lindley_geometrica <- function(xi,theta) # par[1] será theta, par[2] é o p
{  #Restrições: par[1] > 0 e 0 < par[2] < 1
  
  n  <- length(xi)
  f1 <- 2*n*log(theta[1])
  f2 <- n*log(1-theta[2])
  f3 <- sum(log(1+xi))
  f4 <- sum(-theta[1]*xi)
  f5 <- n*log(theta[1]+1)
  f6 <- 2*sum(log(1-theta[2]*(1+theta[1]*xi/(theta[1]+1))*exp(-xi*theta[1])))
  
  ll <- f1 + f2 + f3 + f4 - f5 - f6
  return(ll)
  
}
set.seed(9999)


# 1 tentativa de otimização
x <- rlingley_geom(1000, c(3,0.2))

A <- matrix(c(0, 1, 1,0), 2, 2)           # Restrições dos parâmetros
B <- c(-0.0001, 0.0001)                   # 0theta + 1rho >= 0.0001      e 1theta + 0rho >=0

sann <- maxLik(log_lindley_geometrica,
       start = c(2, 0.3),
       xi=rlingley_geom(1000, c(3,x)),
       method = 'SANN',
       constraints = list(ineqA=A, ineqB=B))

summary(sann)
str(sann)
# O ajuste é horroroso


par_comb <- list(c(1.5, 0.1), c(1.5, 0.5), c(1.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  2, 0.1), c(  2, 0.5), c(  2, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000

simulacoes_sann <- array(c(rep(0,6)), dim=c(N,8,9,10))
simulacoes_sann

A <- matrix(c(0, 1, 1,0), 2, 2)           # Restrições dos parâmetros
B <- c(-0.0001, 0.0001)                   # 0theta + 1rho >= 0.0001      e 1theta + 0rho >=0

help(optim)

set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9)                         # Combinação de parâmetros
  { par <- par_comb[[index_par]]
  amostra <- rlingley_geom(n, par=par)               # Amostra
  op      <- try(maxLik(log_lindley_geometrica,
                        start = c(1.1, 0.3),
                        xi=amostra,
                        method = 'SANN',
                        constraints = list(ineqA=A, ineqB=B)))             
  
  if(typeof(op) == 'character')
  { op      <-  try()                   
  }
  
  if(typeof(op) == 'character')
  { valores <- c(NA, NA, par[1], par[2], n, 99, NA, NA)
  next}
  
  h <- try(solve(op$hessian))              # Tenta inverter a hessiana
  if(typeof(h) == 'character') {h <- c(NA, NA, NA, NA)}  # Se não for invetível, ele guarda o erro em character
  # Daí se o tipo for character, h vira um vetor de NA
  
  valores <- c(op$par[1], op$par[2], par[1], par[2], n, op$convergence, h[1], h[4])
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
  
  cat('itr:', i, '-' , valores, '\n')
  simulacoes_sann[i, ,index_par, index_n] <- valores
  
  }
  }
  
}

simulacoes_sann

# Salvando as simulaçoes

simulacoes_bfgs
simulacoes_cg
simulacoes_nelder
simulacoes_lbfgs
simulacoes_cg[,,8,9]
# save(simulacoes_bfgs, file='simulacoes_bfgs.RData')
# save(simulacoes_cg, file='simulacoes_cg.RData')
# save(simulacoes_nelder, file='simulacoes_nelder.RData')
# save(simulacoes_lbfgs, file='simulacoes_lbfgs.RData')


<<<<<<< HEAD
# Verificando os resultados das simulações --------------------------------


simulacoes_nelder[,,1,2]
simulacoes_nelder[,1,1,] # Todos os theta estimados quando theta = 0.5

simulacoes_nelder[,1,1,][simulacoes_nelder[,1,1,] == 0]

thetas_nelder <- apply(simulacoes_nelder[,1,1,] , MARGIN = 2, FUN= mean) # A média de cada uma das
thetas_nelder
plot(thetas_nelder, type = 'l', ylim = c(0.45, 0.70))
abline(h=0.5, col = 'red')
=======





# Método SANN -------------------------------------------------------------

par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), 
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000

simulacoes_sann <- array(c(rep(0,6)), dim=c(N,8,9,10))







f1 <- function(par, xs) {
  
  n  <- length(xs)
  f1 <- 2*n*log(par[1])
  f2 <- n*log(1-par[2])
  f3 <- sum(log(1+xs))
  f4 <- sum(-par[1]*xs)
  f5 <- n*log(par[1]+1)
  f6 <- 2*sum(log(1-par[2]*(1+par[1]*xs/(par[1]+1))*exp(-xs*par[1])))
  
  ll <- f1 + f2 + f3 + f4 - f5 - f6
  return(ll)
  
  
}


# sa <- function(f, N, rho, t0, theta0) {
#   
#   y <- NULL
#   
#   theta.old <- c()
#   theta.old[1] <- theta0[1]
#   theta.old[2] <- theta0[2]
#   
#   repeat {
#     
#     s <- 0
#     
#     for (k in 1:N) {
#       f.old <- f(theta.old, xs)
#       theta.new <- runif(1, -1, 1) + t
#     }
#     
#   }
#   
# }


# GABRIEL FAZENDO
library(maxLik)
?maxLik
for (i in 1:N) {
  for (index_n in 1:10) {
    
    
    n <- seq(1, 100, 10)[index_n]
    
    for (index_par in 1:9) {
      
      par <- par_comb[[index_par]]
      amostra <- rlingley_geom(n, par=par)
      op <- try(maxLik(logLik = log_lindley_geometrica, start = c(1.5, 0.3), 
                       method = 'SANN', x = amostra))
      
      
      
    }
    
    h <- try(solve(op$hessian))              # Tenta inverter a hessiana
    if(typeof(h) == 'character') {h <- c(NA, NA, NA, NA)}
    
    valores <- c(coef(op)[1], coef(op)[2], par[1], par[2], n, 1,  h[1], h[4])
    # Valores recebe o que queremos dessa bagaça toda,
    # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
    
    cat('itr:', i, '-' , valores, '\n')
    simulacoes_lbfgs[i, ,index_par, index_n] <- valores
  }
  
  
}

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
>>>>>>> eb4dc0105d987a8ae7aebdd973062335bae5eebb
