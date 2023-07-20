
# Optim -------------------------------------------------------------------
#EXEMPLO
x <- rexp(100, rate=5)

emv <- 1/mean(x)
emv

hist(x, freq = F)
curve(dexp(x,rate = emv), add = T,lwd = 2)

x <- rexp(1000, rate =,)
emv <- 1/mean(x)
hist(X, freq = F)
curve(dexp(x, rate = ))

log.like.exp <- function(lambda, data)
{
  n <- length(data) 
  ll <- n*log(lambda) - lambda * n * mean(data)
  return(ll)
}
x <- rexp(1000, rate =2)
emv.analitico <- 1/mean(x)
emv.numerico <- optim(par = 2, fn = log.like.exp, data = x, control = list(fnscale = -1))

emv.analitico
emv.numerico   


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



log_lindley_geometrica <- function(x, n, par) # par[1] será theta, par[2] é o p
{
  2*n*log(par[1]) + n*log(1-par[2]) + sum(log(1+x)) + 
    sum(log(exp(-theta[1]*x))) - n*log(par[1] + 1) -
    2*sum(log(1-theta[2]*(1+(theta[1]*x))))
  
  
  f1 <- 2*n*log(par[1])
  f2 <- n*log(1-theta[2])
  f3 <- sum(log(1+x))
  f4 <- sum(log(exp(-theta[1]*x)))
  f5 <- n*log(theta[1]+1)
  f6 <- 2*sum(log(1-theta[2]*(1+theta[1]*x/(theta[1]+1))*exp(-x*par[1])))
  
  ll <- f1 + f2 + f3 + f4 - f5 - f6
  return(ll)
  
}


