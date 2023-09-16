library(LambertW)

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

lindley_geom <- function(u, par)
  
{
  e           <- exp(-par[1] - 1)
  p1          <- ((u - 1) * (par[1] + 1) * e) / (1 - par[2]*u)
  quantile    <- -1 - 1 / par[1] - W(p1, -1)* par[1]^(-1)
  
  return(quantile)
}


rlindley_geom_AA <- function(n, par)
{
  
  u <- runif(n/2, min = 0, max = 1) 
  values     <- qlindley_geom(u, par) 
  values_inv <- qlindley_geom(1 - u, par) 
  return(list(amostra1 = values, amostra2 = values_inv))
  
}



# Método SANN -----------------------------------------------------------

library(LambertW)


par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50

AA_simulacoes_sann <- array(c(rep(0,6)), dim=c(N,13,9,10)) # COLUNAS, theta_op1, theta_op2, rho_op1, rho_op2, par1, par2, n, var_theta_op1, var_theta_op2, var_rho_op1, var_rho_op2, cov_theta, cov_rho
AA_simulacoes_sann

set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{ 
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9)                         # Combinação de parâmetros
  {
  par <- par_comb[[index_par]]
  
  correlacao <- 1
  while(correlacao > 0)
  { amostra <- rlindley_geom_AA(n, par=par)   # Amostra
  amostra1 <- amostra$amostra1
  amostra2 <- amostra$amostra2
  correlacao <- cor(amostra1, amostra2)}
  op1      <-  try(optim(par = par,            # Chute inicial
                         fn = log_lindley_geometrica,  # Log-Verossimilhança
                         x = amostra1,                  # Amostra
                         control = list(fnscale = -1),
                         method = 'SANN'))                   
  
  
  op2      <-  try(optim(par = par,            # Chute inicial
                         fn = log_lindley_geometrica,  # Log-Verossimilhança
                         x = amostra2,                  # Amostra
                         control = list(fnscale = -1),
                         method = 'SANN'))                   
  
  
  if(typeof(op2) == 'character' || typeof(op1) == 'character')
  { valores <- c(NA, NA, NA, NA, par[1], par[2], n, rep(NA,6))
  AA_simulacoes_sann[i, ,index_par, index_n] <- valores
  next}
  
  
  valores <- c(op1$par[1], op2$par[1],op1$par[2], op2$par[2], par[1], par[2], n, rep(NA,6))
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
  
  cat('itr:', i, '-' , valores, '\n')
  AA_simulacoes_sann[i, ,index_par, index_n] <- valores
  
  }
  }
  
}

for(i in 1:10)
{
  for(j in 1:9)
  {
    matriz_theta <- cov(na.omit(AA_simulacoes_sann[,c(1,2),j,i]))
    matriz_rho   <- cov(na.omit(AA_simulacoes_sann[,c(3,4),j,i]))
    
    
    AA_simulacoes_sann[,8,j,i]  <- matriz_theta[1,1]
    AA_simulacoes_sann[,9,j,i]  <- matriz_theta[2,2]
    AA_simulacoes_sann[,10,j,i] <- matriz_rho[1,1]
    AA_simulacoes_sann[,11,j,i] <- matriz_rho[2,2]
    
    AA_simulacoes_sann[,12,j,i] <- matriz_theta[1,2]
    AA_simulacoes_sann[,13,j,i]<- matriz_rho[1,2]
  }
}

AA_simulacoes_sann
save(AA_simulacoes_sann, file= 'AA_simulacoes_sann.RData')
AA_dignostico_sann <- array(c(rep(0,6)), dim=c(N,13,9,10))   
# n, Var_theta, Var_rho, vicio_theta, vicio_rho, eqm_theta, eqm_rho, erro_theta, erro_rho, prob_theta, prob_rho, tam_theta, tam_rho,

AA_dignostico_sann[,1,,]  <- AA_simulacoes_sann[,7,,] # n
AA_dignostico_sann[,2,,]  <- (AA_simulacoes_sann[,8,,] + AA_simulacoes_sann[,9,,])*0.25 + 0.5*AA_simulacoes_sann[,c(12),,] # var theta
AA_dignostico_sann[,3,,]  <-  (AA_simulacoes_sann[,10,,] + AA_simulacoes_sann[,11,,])*0.25 + 0.5*AA_dignostico_sann[,c(13),,] # var rho
AA_dignostico_sann[,4,,]  <- (AA_simulacoes_sann[,1,,] + AA_simulacoes_sann[,2,,])*0.5 - AA_simulacoes_sann[,5,,]     # Vício de theta
AA_dignostico_sann[,5,,]  <- (AA_simulacoes_sann[,3,,] + AA_simulacoes_sann[,4,,])*0.5 - AA_simulacoes_sann[,6,,]     # Vício de rho
AA_dignostico_sann[,6,,]  <- ((AA_simulacoes_sann[,1,,] + AA_simulacoes_sann[,2,,])*0.5 - AA_simulacoes_sann[,5,,])^2 # EQM de theta
AA_dignostico_sann[,7,,]  <- ((AA_simulacoes_sann[,3,,] + AA_simulacoes_sann[,4,,])*0.5 - AA_simulacoes_sann[,6,,])^2 # EQM de rho
AA_dignostico_sann[,8,,]  <- sqrt(AA_simulacoes_sann[,2,,])  # Erro theta
AA_dignostico_sann[,9,,]  <- sqrt(AA_simulacoes_sann[,3,,])  # Erro rho

AA_dignostico_sann[,10,,] <- (( mean(AA_simulacoes_sann[,c(1,2),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_sann[,8,,]) < AA_simulacoes_sann[,5,,]) & (( mean(AA_simulacoes_sann[,c(1,2),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_sann[,8,,]) > AA_simulacoes_sann[,5,,]) # Prob cobertura theta

AA_dignostico_sann[,11,,] <- (( mean(AA_simulacoes_sann[,c(3,4),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_sann[,9,,]) < AA_simulacoes_sann[,6,,]) & (( mean(AA_simulacoes_sann[,c(3,4),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_sann[,9,,]) > AA_simulacoes_sann[,6,,]) # Prob cobertura rho


AA_dignostico_sann[,12,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_sann[,8,,] # Tam_t
AA_dignostico_sann[,13,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_sann[,9,,] # Tam_r
AA_dignostico_sann



save(AA_simulacoes_sann,  file = 'AA_simulacoes_sann.Rdata')
save(AA_dignostico_sann, file = 'AA_dignostico_sann.Rdata')
