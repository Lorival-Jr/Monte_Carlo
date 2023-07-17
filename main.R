# Geração de Valores Aleatórios -------------------------------------------


# 9. Geração de Valores Aleatórios

# ENTENDENDO O CÓDIGO
my.dexp             <-      function(x, lambda)
{
  lambda*exp(-lambda*x) # Criando a função densidade
}

my.pexp.numerical   <-      function(x, lambda)
{
  my.int          <-      function(x, lambda) # Criando uma função pra utilizar no apply
  { # Ela integra a fdp, recebendo o parametro lambda pra fdp e o lower sendo o intervalo de integração inf = 0 e o sup = x
    integrate(my.dexp, lambda = lambda, lower = 0, upper = x)$value # é pego do resultado apenas o valor
  } 
  sapply(x, FUN=my.int, lambda) # Aplica nos dados a integral da fdp, ou seja a fda
}

my.qexp.numerical   <-      function(q, lambda) # Função quantil
{ 
  f                 <-      function(P, fixed)  # Função pra ser aplicada no unirrot                    
  {
    lambda          <-      fixed$lambda        # fixed é uma lista com 2 dimensões
    q               <-      fixed$q             # A 1° o valor de lambda, e a 2° todos quantis que queremos calcular
    criterion       <-      q - my.pexp.numerical(P, lambda)   # É o inverso da probabiladade da FDA                     
    return(criterion)
  }
  P                 <-      numeric(length(q)) # Para cada q passado vai achar uma raiz                                     
  for(i in 1:length(q))
  {
    fixed           <-      list(lambda = lambda, q = q[i]) # É a lista que será passado para função f
    root.p          <-      uniroot(f, lower = 0, upper = 100, fixed = fixed) 
    # unirrot procura as raizes da função f, com o intervalo de integração de 0 a 100, o fixed é a lista que a f recebe
    P[i]            <-      root.p$root # pega apenas a raiz pra retornar
  }
  return(P)
}

my.qexp.numerical(0.9, 2)
my.qexp.numerical(c(0.9, 0.1), 2)


my.rexp <- function(n, lambda)
{
  u <- runif(n, min = 0, max = 1) # pega um valor aleatorio entre 0 e 1 e repassa ele pra funcao quantil criada
  rnd.values <- my.qexp(u, lambda) # A função quantil usa a probabilidade gerada e retorna um quantil
  return(rnd.values)
}

# Lindley Geométrica ------------------------------------------------------

### Forma analítica

dlindley_geom <- function(x, theta, p)
{
  e           <- exp(-theta * x) 
  p1          <- (theta^2 / (theta + 1)) * (1-p) * (1+ x) * e
  p2          <- (1 - p*(1 + (theta * x)/(theta + 1)) * e)^(-2)
  pdf         <- p1 * p2
  
  return(pdf)
}
# restrições theta > 0 e 0 < p < 1
#  stopifnot(p < 1, p > 0, theta > 0, x > 0) no R

plindley_geom <- function(x, theta, p) 
{
  e              <- exp(-theta * x)
  p1             <- (1 + (theta*x)/(theta + 1))*e
  cdf            <- (1 - p1)/(1 - p*p1) 
  
  return(cdf)
}

#install.packages('LambertW')
library('LambertW')
??LambertW

qlindley_geom <- function(u, theta, p)
  
{
  e           <- exp(-theta - 1)
  p1          <- ((u - 1) * (theta + 1) * e) / (1 - p*u)
  quantile    <- -1 - 1 / theta - W(p1, -1)* theta^(-1)
  
  return(quantile)
}

rlingley_geom <- function(n, theta, p)
{
  
  u <- runif(n, min = 0, max = 1) 
  rnd.values <- qlindley_geom(u, theta, p) 
  return(rnd.values)
  
}
## Testando

#install.packages("LindleyPowerSeries")
library("LindleyPowerSeries")

# Pagina 4
# https://cran.r-project.org/web/packages/LindleyPowerSeries/LindleyPowerSeries.pdf
help(dlindleygeometric)

dlindley_geom(c(1, 1.5, 4), 1, 0.4)
dlindleygeometric(c(1, 1.5, 4), 2, 0.4) # lambda = theta e theta = p


plindley_geom(c(1, 1.5, 4), 1, 0.4)
plindleygeometric(c(1,1.5,4), 2, 0.4)


qlindley_geom(c(0.3,0.5,0.8), 2, 0.4)
qlindleygeometric(c(0.3,0.5,0.8), 2, 0.4)


rlingley_geom(100, 2, 0.4)
rlindleygeometric(100, 2, 0.4)


par(mfrow = c(1,2))

hist(rlingley_geom(1000, 0.8, 0.4), freq = F)
hist(rlindleygeometric(1000, 0.8, 0.4), freq = F)

hist(rlingley_geom(1000, 0.2, 0.4), freq = F)
hist(rgeom(1000, 0.4), freq = F)
par(mfrow = c(1,1))

hist(rlingley_geom(1000, 0.8, 0.4), freq = F)
curve(dlindleygeometric(x, 0.8, 0.4), add=T, col = 'green', xlim= c(0.001,13), lwd = 2)

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
