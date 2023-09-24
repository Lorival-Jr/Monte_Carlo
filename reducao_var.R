# Redução de Variância ----------------------------------------------------

# A redução se dará via amostragem antitética

library(LambertW)

rlindley_geom_AA <- function(n, par)
{
  
  u <- runif(n/2, min = 0, max = 1) 
  values     <- qlindley_geom(u, par) 
  values_inv <- qlindley_geom(1 - u, par) 
  return(list(amostra1 = values, amostra2 = values_inv))
  
}
cov(amostra[[1]], amostra[[2]])
amostra <- rlindley_geom_AA(20, c(1,0.5));cor(amostra[[1]], amostra[[2]])

log_lindley_geometrica2 <- function(xi, par) # par[1] será theta, par[2] é o p
{  #Restrições: par[1] > 0 e 0 < par[2] < 1
  
  
  # 2*n*log(par[1]) + n*log(1-par[2]) + sum(log(1+x)) + 
  #   sum(log(exp(-theta[1]*x))) - n*log(par[1] + 1) -
  #   2*sum(log(1-theta[2]*(1+(theta[1]*x))))
  
  n  <- length(xi)
  f1 <- 2*n*log(par[1])
  f2 <- n*log(1-par[2])
  f3 <- sum(log(1+xi))
  f4 <- sum(-par[1]*xi)
  f5 <- n*log(par[1]+1)
  f6 <- 2*sum(log(1-par[2]*(1+par[1]*xi/(par[1]+1))*exp(-xi*par[1])))
  
  ll <- f1 + f2 + f3 + f4 - f5 - f6
  return(ll)
  
}



maxLik(log_lindley_geometrica2,
       start = c(0.5,0.5),
       xi=amostra,
       method = 'NR')

par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8),   # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),   # Daí são as combinações
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000                                                # O N de simulações pedido é 50000

AA_simulacoes_newton <- array(c(rep(0,6)), dim=c(N,13,9,10))  # Esse array vai guardar os resultados
AA_simulacoes_newton                                         # Basicamente são 90 matrizes dentro de um array

dim(AA_simulacoes_newton) # Serão 50 mil linhas, 15 colunas e 90 matrizes



library(maxLik)
set.seed(9999)
for (i in 1:N) # Número de simulações
{ 
  for (index_n in 1:10) # Tamanho da amostra
  {n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9) # Combinação de parâmetros
  { par <- par_comb[[index_par]]
 
  
  correlacao <- 1
  while(correlacao > 0)
  { amostra <- rlindley_geom_AA(n, par=par)   # Amostra
  amostra1 <- amostra$amostra1
  amostra2 <- amostra$amostra2
  correlacao <- cor(amostra1, amostra2)}
  
  op1 <-       try(maxLik(log_lindley_geometrica2,
                         start = c(par),
                         xi=amostra1,
                         method = 'NR'))
  
  op2 <-       try(maxLik(log_lindley_geometrica2,
                          start = c(par),
                          xi=amostra2,
                          method = 'NR'))

  
  if(typeof(op2) == 'character' || typeof(op1) == 'character')
  { valores <- c(NA, NA, NA, NA, par[1], par[2], n, rep(NA,6))
  AA_simulacoes_newton[i, ,index_par, index_n] <- valores
  next}


  valores <- c(op1$estimate[1], op2$estimate[1], op1$estimate[2], op2$estimate[2], par[1], par[2], n, rep(NA,6))

  
  cat('itr:', i, '-' , valores, '\n')    # Inútil, é só pra vc ficar vendo oq ta acontecendo
  
  AA_simulacoes_newton[i, ,index_par, index_n] <- valores  # Guarda na tabela
  
  }
  }
  
}

for(i in 1:10)
{
  for(j in 1:9)
  {
    matriz_theta <- cov(na.omit(AA_simulacoes_newton[,c(1,2),j,i]))
    matriz_rho   <- cov(na.omit(AA_simulacoes_newton[,c(3,4),j,i]))
    
    
    AA_simulacoes_newton[,8,j,i]  <- matriz_theta[1,1]
    AA_simulacoes_newton[,9,j,i]  <- matriz_theta[2,2]
    AA_simulacoes_newton[,10,j,i] <- matriz_rho[1,1]
    AA_simulacoes_newton[,11,j,i] <- matriz_rho[2,2]
    
    AA_simulacoes_newton[,12,j,i] <- matriz_theta[1,2]
    AA_simulacoes_newton[,13,j,i]<- matriz_rho[1,2]
  }
}





save(AA_simulacoes_newton, file='AA_simulacoes_newton.RData')

load('AA_simulacoes_newton.RData')

AA_dignostico_newton <- array(c(rep(0,6)), dim=c(N,13,9,10))   
# n, Var_theta, Var_rho, vicio_theta, vicio_rho, eqm_theta, eqm_rho, erro_theta, erro_rho, prob_theta, prob_rho, tam_theta, tam_rho,

AA_dignostico_newton[,1,,]  <- AA_simulacoes_newton[,7,,]
AA_dignostico_newton[,2,,]  <- (AA_simulacoes_newton[,8,,] + AA_simulacoes_newton[,9,,])*0.25 + 0.5*AA_simulacoes_newton[,c(12),,]
AA_dignostico_newton[,3,,]  <-  (AA_simulacoes_newton[,10,,] + AA_simulacoes_newton[,11,,])*0.25 + 0.5*AA_simulacoes_newton[,c(13),,]
AA_dignostico_newton[,4,,]  <- (AA_simulacoes_newton[,1,,] + AA_simulacoes_newton[,2,,])*0.5 - AA_simulacoes_newton[,5,,]     # Vício de theta
AA_dignostico_newton[,5,,]  <- (AA_simulacoes_newton[,3,,] + AA_simulacoes_newton[,4,,])*0.5 - AA_simulacoes_newton[,6,,]     # Vício de rho
AA_dignostico_newton[,6,,]  <- ((AA_simulacoes_newton[,1,,] + AA_simulacoes_newton[,2,,])*0.5 - AA_simulacoes_newton[,5,,])^2 # EQM de rho
AA_dignostico_newton[,7,,]  <- ((AA_simulacoes_newton[,3,,] + AA_simulacoes_newton[,4,,])*0.5 - AA_simulacoes_newton[,6,,])^2 # EQM de rho
AA_dignostico_newton[,8,,]  <- sqrt(AA_simulacoes_newton[,2,,])  # Erro theta
AA_dignostico_newton[,9,,]  <- sqrt(AA_simulacoes_newton[,3,,])  # Erro rho
AA_dignostico_newton[,10,,] <- (( mean(AA_simulacoes_newton[,c(1,2),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_newton[,8,,]) < AA_simulacoes_newton[,5,,]) & (( mean(AA_simulacoes_newton[,c(1,2),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_newton[,8,,]) > AA_simulacoes_newton[,5,,]) # Prob cobertura theta
AA_dignostico_newton[,11,,] <- (( mean(AA_simulacoes_newton[,c(3,4),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_newton[,9,,]) < AA_simulacoes_newton[,6,,]) & (( mean(AA_simulacoes_newton[,c(3,4),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_newton[,9,,]) > AA_simulacoes_newton[,6,,]) # Prob cobertura rho
AA_dignostico_newton[,12,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_newton[,8,,] # Tam_t
AA_dignostico_newton[,13,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_newton[,9,,] # Tam_r
AA_dignostico_newton

save(AA_dignostico_newton, file='AA_dignostico_newton.RData')



# Método NELDER ---------------------------------------------------------

library(LambertW)


par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000

AA_simulacoes_nelder <- array(c(rep(0,6)), dim=c(N,13,9,10)) # COLUNAS, theta_op1, theta_op2, rho_op1, rho_op2, par1, par2, n, var_theta_op1, var_theta_op2, var_rho_op1, var_rho_op2, cov_theta, cov_rho
AA_simulacoes_nelder


set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{ 
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9)                         # Combinação de parâmetros
  { par <- par_comb[[index_par]]
  
  correlacao <- 1
  while(correlacao > 0)
  { amostra <- rlindley_geom_AA(n, par=par)   # Amostra
  amostra1 <- amostra$amostra1
  amostra2 <- amostra$amostra2
  correlacao <- cor(amostra1, amostra2)}
  
  
  
  op1      <- try(optim(par = par,            # Chute inicial
                        fn = log_lindley_geometrica,  # Log-Verossimilhança
                        x = amostra1,                  # Amostra
                        control = list(fnscale = -1),
                        method = 'Nelder-Mead',          # Método
                        hessian = F
  ))                 
  
  op2      <- try(optim(par = par,            # Chute inicial
                        fn = log_lindley_geometrica,  # Log-Verossimilhança
                        x = amostra2,                  # Amostra
                        control = list(fnscale = -1),
                        method = 'Nelder-Mead',          # Método
                        hessian = F                 # Calcular a hessiana
                        
  ))                 
  
  
  if(typeof(op2) == 'character' || typeof(op1) == 'character')
  { valores <- c(NA, NA, NA, NA, par[1], par[2], n, rep(NA,6))
  AA_simulacoes_nelder[i, ,index_par, index_n] <- valores
  next}
  
  
  
  valores <- c(op1$par[1], op2$par[1],op1$par[2], op2$par[2], par[1], par[2], n, rep(NA,6))
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
  
  cat('itr:', i, '-' , valores, '\n')
  AA_simulacoes_nelder[i, ,index_par, index_n] <- valores
  
  }
  }
  
}
save(AA_simulacoes_nelder, file = 'AA_simulacoes_nelder.RData')
load('AA_simulacoes_nelder.RData')
for(i in 1:10)
{
  for(j in 1:9)
  {
    matriz_theta <- cov(na.omit(AA_simulacoes_nelder[,c(1,2),j,i]))
    matriz_rho   <- cov(na.omit(AA_simulacoes_nelder[,c(3,4),j,i]))
    
    
    AA_simulacoes_nelder[,8,j,i]  <- matriz_theta[1,1]
    AA_simulacoes_nelder[,9,j,i]  <- matriz_theta[2,2]
    AA_simulacoes_nelder[,10,j,i] <- matriz_rho[1,1]
    AA_simulacoes_nelder[,11,j,i] <- matriz_rho[2,2]
    
    AA_simulacoes_nelder[,12,j,i] <- matriz_theta[1,2]
    AA_simulacoes_nelder[,13,j,i]<- matriz_rho[1,2]
  }
}



AA_simulacoes_nelder

AA_dignostico_nelder <- array(c(rep(0,6)), dim=c(N,13,9,10))   
# n, Var_theta, Var_rho, vicio_theta, vicio_rho, eqm_theta, eqm_rho, erro_theta, erro_rho, prob_theta, prob_rho, tam_theta, tam_rho,

AA_dignostico_nelder[,1,,]  <- AA_simulacoes_nelder[,7,,]
AA_dignostico_nelder[,2,,]  <- (AA_simulacoes_nelder[,8,,] + AA_simulacoes_nelder[,9,,])*0.25 + 0.5*AA_simulacoes_nelder[,c(12),,]
AA_dignostico_nelder[,3,,]  <-  (AA_simulacoes_nelder[,10,,] + AA_simulacoes_nelder[,11,,])*0.25 + 0.5*AA_simulacoes_nelder[,c(13),,]
AA_dignostico_nelder[,4,,]  <- (AA_simulacoes_nelder[,1,,] + AA_simulacoes_nelder[,2,,])*0.5 - AA_simulacoes_nelder[,5,,]     # Vício de theta
AA_dignostico_nelder[,5,,]  <- (AA_simulacoes_nelder[,3,,] + AA_simulacoes_nelder[,4,,])*0.5 - AA_simulacoes_nelder[,6,,]     # Vício de rho
AA_dignostico_nelder[,6,,]  <- ((AA_simulacoes_nelder[,1,,] + AA_simulacoes_nelder[,2,,])*0.5 - AA_simulacoes_nelder[,5,,])^2 # EQM de rho
AA_dignostico_nelder[,7,,]  <- ((AA_simulacoes_nelder[,3,,] + AA_simulacoes_nelder[,4,,])*0.5 - AA_simulacoes_nelder[,6,,])^2 # EQM de rho
AA_dignostico_nelder[,8,,]  <- sqrt(AA_simulacoes_nelder[,2,,])  # Erro theta
AA_dignostico_nelder[,9,,]  <- sqrt(AA_simulacoes_nelder[,3,,])  # Erro rho
AA_dignostico_nelder[,10,,] <- (( mean(AA_simulacoes_nelder[,c(1,2),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_nelder[,8,,]) < AA_simulacoes_nelder[,5,,]) & (( mean(AA_simulacoes_nelder[,c(1,2),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_nelder[,8,,]) > AA_simulacoes_nelder[,5,,]) # Prob cobertura theta
AA_dignostico_nelder[,11,,] <- (( mean(AA_simulacoes_nelder[,c(3,4),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_nelder[,9,,]) < AA_simulacoes_nelder[,6,,]) & (( mean(AA_simulacoes_nelder[,c(3,4),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_nelder[,9,,]) > AA_simulacoes_nelder[,6,,]) # Prob cobertura rho
AA_dignostico_nelder[,12,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_nelder[,8,,] # Tam_t
AA_dignostico_nelder[,13,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_nelder[,9,,] # Tam_r
AA_dignostico_nelder



save(AA_simulacoes_nelder,  file = 'AA_simulacoes_nelder.Rdata')
save(AA_dignostico_nelder, file = 'AA_dignostico_nelder.Rdata')
load('AA_simulacoes_nelder.Rdata')






# Método BFGS ---------------------------------------------------------

library(LambertW)


par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000

AA_simulacoes_bfgs <- array(c(rep(0,6)), dim=c(N,13,9,10)) # COLUNAS, theta_op1, theta_op2, rho_op1, rho_op2, par1, par2, n, var_theta_op1, var_theta_op2, var_rho_op1, var_rho_op2, cov_theta, cov_rho
AA_simulacoes_bfgs


set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{ 
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9)                         # Combinação de parâmetros
  { par <- par_comb[[index_par]]
  
  correlacao <- 1
  while(correlacao > 0)
  { amostra <- rlindley_geom_AA(n, par=par)   # Amostra
  amostra1 <- amostra$amostra1
  amostra2 <- amostra$amostra2
  correlacao <- cor(amostra1, amostra2)}
  
  
  
  op1      <- try(optim(par = par,            # Chute inicial
                        fn = log_lindley_geometrica,  # Log-Verossimilhança
                        x = amostra1,                  # Amostra
                        control = list(fnscale = -1),
                        method = 'BFGS',              # Método
                        hessian = F                 # Calcular a hessiana
  ))                 
  
  op2      <- try(optim(par = par,            # Chute inicial
                        fn = log_lindley_geometrica,  # Log-Verossimilhança
                        x = amostra2,                  # Amostra
                        control = list(fnscale = -1),
                        method = 'BFGS',          # Método
                        hessian = F                # Calcular a hessiana
                        
  ))                 
  
  if(typeof(op2) == 'character' || typeof(op1) == 'character')
  { valores <- c(NA, NA, NA, NA, par[1], par[2], n, rep(NA,6))
  AA_simulacoes_lbfgs[i, ,index_par, index_n] <- valores
  next}
  
  valores <- c(op1$par[1], op2$par[1],op1$par[2], op2$par[2], par[1], par[2], n, rep(NA,6))
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
  
  cat('itr:', i, '-' , valores, '\n')
  AA_simulacoes_bfgs[i, ,index_par, index_n] <- valores
  
  }
  }
  
}
for(i in 1:10)
{
  for(j in 1:9)
  {
    matriz_theta <- cov(na.omit(AA_simulacoes_bfgs[,c(1,2),j,i]))
    matriz_rho   <- cov(na.omit(AA_simulacoes_bfgs[,c(3,4),j,i]))
    
    
    AA_simulacoes_bfgs[,8,j,i]  <- matriz_theta[1,1]
    AA_simulacoes_bfgs[,9,j,i]  <- matriz_theta[2,2]
    AA_simulacoes_bfgs[,10,j,i] <- matriz_rho[1,1]
    AA_simulacoes_bfgs[,11,j,i] <- matriz_rho[2,2]
    
    AA_simulacoes_bfgs[,12,j,i] <- matriz_theta[1,2]
    AA_simulacoes_bfgs[,13,j,i]<- matriz_rho[1,2]
  }
}



AA_simulacoes_bfgs
save(AA_simulacoes_bfgs, file = 'RDatas/AA_simulacoes_bfgs.RData')
AA_dignostico_bfgs <- array(c(rep(0,6)), dim=c(N,13,9,10))   
# n, Var_theta, Var_rho, vicio_theta, vicio_rho, eqm_theta, eqm_rho, erro_theta, erro_rho, prob_theta, prob_rho, tam_theta, tam_rho,
load('RDatas/AA_simulacoes_bfgs.RData')
AA_dignostico_bfgs[,1,,]  <- AA_simulacoes_lbfgs[,7,,]
AA_dignostico_bfgs[,2,,]  <- (AA_simulacoes_lbfgs[,8,,] + AA_simulacoes_lbfgs[,9,,])*0.25 + 0.5*AA_simulacoes_lbfgs[,c(12),,]
AA_dignostico_bfgs[,3,,]  <-  (AA_simulacoes_lbfgs[,10,,] + AA_simulacoes_lbfgs[,11,,])*0.25 + 0.5*AA_simulacoes_lbfgs[,c(13),,]
AA_dignostico_bfgs[,4,,]  <- (AA_simulacoes_lbfgs[,1,,] + AA_simulacoes_lbfgs[,2,,])*0.5 - AA_simulacoes_lbfgs[,5,,]     # Vício de theta
AA_dignostico_bfgs[,5,,]  <- (AA_simulacoes_lbfgs[,3,,] + AA_simulacoes_lbfgs[,4,,])*0.5 - AA_simulacoes_lbfgs[,6,,]     # Vício de rho
AA_dignostico_bfgs[,6,,]  <- ((AA_simulacoes_lbfgs[,1,,] + AA_simulacoes_lbfgs[,2,,])*0.5 - AA_simulacoes_lbfgs[,5,,])^2 # EQM de rho
AA_dignostico_bfgs[,7,,]  <- ((AA_simulacoes_lbfgs[,3,,] + AA_simulacoes_lbfgs[,4,,])*0.5 - AA_simulacoes_lbfgs[,6,,])^2 # EQM de rho
AA_dignostico_bfgs[,8,,]  <- sqrt(AA_simulacoes_lbfgs[,2,,])  # Erro theta
AA_dignostico_bfgs[,9,,]  <- sqrt(AA_simulacoes_lbfgs[,3,,])  # Erro rho
AA_dignostico_bfgs[,10,,] <- (( mean(AA_simulacoes_lbfgs[,c(1,2),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_bfgs[,8,,]) < AA_simulacoes_lbfgs[,5,,]) & (( mean(AA_simulacoes_lbfgs[,c(1,2),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_bfgs[,8,,]) > AA_simulacoes_lbfgs[,5,,]) # Prob cobertura theta
AA_dignostico_bfgs[,11,,] <- (( mean(AA_simulacoes_lbfgs[,c(3,4),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_bfgs[,9,,]) < AA_simulacoes_lbfgs[,6,,]) & (( mean(AA_simulacoes_lbfgs[,c(3,4),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_bfgs[,9,,]) > AA_simulacoes_lbfgs[,6,,]) # Prob cobertura rho
AA_dignostico_bfgs[,12,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_bfgs[,8,,] # Tam_t
AA_dignostico_bfgs[,13,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_bfgs[,9,,] # Tam_r
AA_dignostico_bfgs



save(AA_simulacoes_bfgs,  file = 'AA_simulacoes_bfgs.Rdata')
save(AA_dignostico_bfgs, file = 'AA_dignostico_bfgs.Rdata')

AA_simulacoes_lbfgs[,7,1,]
simulacoes_bfgs[,7,,]




AA_simulacoes_lbfgs[,,2,3]

# Método L-BFGS-B ---------------------------------------------------------

library(LambertW)


par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50

AA_simulacoes_lbfgs <- array(c(rep(0,6)), dim=c(N,13,9,10)) # COLUNAS, theta_op1, theta_op2, rho_op1, rho_op2, par1, par2, n, var_theta_op1, var_theta_op2, var_rho_op1, var_rho_op2, cov_theta, cov_rho
AA_simulacoes_lbfgs


set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{ 
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9)                         # Combinação de parâmetros
  { par <- par_comb[[index_par]]
  amostra <- rlindley_geom_AA(n, par=par)   # Amostra
  amostra1 <- amostra$amostra1
  amostra2 <- amostra$amostra2
  
  op1      <-  try(optim(par = par,            # Chute inicial
                         fn = log_lindley_geometrica,  # Log-Verossimilhança
                         x = amostra1,                  # Amostra
                         control = list(fnscale = -1),
                         method = 'L-BFGS-B',          # Método
                         lower = c(0.000001,  0.000001),
                         upper = c(100, 0.999999)
  ))                   
  
  
  op2      <-  try(optim(par = par,            # Chute inicial
                         fn = log_lindley_geometrica,  # Log-Verossimilhança
                         x = amostra2,                  # Amostra
                         control = list(fnscale = -1),
                         method = 'L-BFGS-B',          # Método
                         lower = c(0.000001,  0.000001),
                         upper = c(100, 0.999999)
  ))                   
  
  
  if(typeof(op2) == 'character' || typeof(op1) == 'character')
  { valores <- c(NA, NA, NA, NA, par[1], par[2], n, rep(NA,6))
  AA_simulacoes_lbfgs[i, ,index_par, index_n] <- valores
  next}
  
  
  valores <- c(op1$par[1], op2$par[1],op1$par[2], op2$par[2], par[1], par[2], n, rep(NA,6))
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
  
  cat('itr:', i, '-' , valores, '\n')
  AA_simulacoes_lbfgs[i, ,index_par, index_n] <- valores
  
  }
  }
  
}
load('AA_simulacoes_lbfgs.Rdata')
cov(na.omit(AA_simulacoes_lbfgs[,c(1,2),1,1]))

for(i in 1:10)
{
  for(j in 1:9)
  {
    matriz_theta <- cov(na.omit(AA_simulacoes_lbfgs[,c(1,2),j,i]))
    matriz_rho   <- cov(na.omit(AA_simulacoes_lbfgs[,c(3,4),j,i]))
    
    
    AA_simulacoes_lbfgs[,8,j,i]  <- matriz_theta[1,1]
    AA_simulacoes_lbfgs[,9,j,i]  <- matriz_theta[2,2]
    AA_simulacoes_lbfgs[,10,j,i] <- matriz_rho[1,1]
    AA_simulacoes_lbfgs[,11,j,i] <- matriz_rho[2,2]
    
    AA_simulacoes_lbfgs[,12,j,i] <- matriz_theta[1,2]
    AA_simulacoes_lbfgs[,13,j,i]<- matriz_rho[1,2]
  }
}

AA_simulacoes_lbfgs
save(AA_simulacoes_lbfgs, file= 'AA_simulacoes_lbfgs.RData')
AA_dignostico_lbfgs <- array(c(rep(0,6)), dim=c(N,13,9,10))   
# n, Var_theta, Var_rho, vicio_theta, vicio_rho, eqm_theta, eqm_rho, erro_theta, erro_rho, prob_theta, prob_rho, tam_theta, tam_rho,

AA_dignostico_lbfgs[,1,,]  <- AA_simulacoes_lbfgs[,7,,] # n
AA_dignostico_lbfgs[,2,,]  <- (AA_simulacoes_lbfgs[,8,,] + AA_simulacoes_lbfgs[,9,,])*0.25 + 0.5*AA_simulacoes_lbfgs[,c(12),,] # var theta
AA_dignostico_lbfgs[,3,,]  <-  (AA_simulacoes_lbfgs[,10,,] + AA_simulacoes_lbfgs[,11,,])*0.25 + 0.5*AA_simulacoes_lbfgs[,c(13),,] # var rho
AA_dignostico_lbfgs[,4,,]  <- (AA_simulacoes_lbfgs[,1,,] + AA_simulacoes_lbfgs[,2,,])*0.5 - AA_simulacoes_lbfgs[,5,,]     # Vício de theta
AA_dignostico_lbfgs[,5,,]  <- (AA_simulacoes_lbfgs[,3,,] + AA_simulacoes_lbfgs[,4,,])*0.5 - AA_simulacoes_lbfgs[,6,,]     # Vício de rho
AA_dignostico_lbfgs[,6,,]  <- ((AA_simulacoes_lbfgs[,1,,] + AA_simulacoes_lbfgs[,2,,])*0.5 - AA_simulacoes_lbfgs[,5,,])^2 # EQM de theta
AA_dignostico_lbfgs[,7,,]  <- ((AA_simulacoes_lbfgs[,3,,] + AA_simulacoes_lbfgs[,4,,])*0.5 - AA_simulacoes_lbfgs[,6,,])^2 # EQM de rho
AA_dignostico_lbfgs[,8,,]  <- sqrt(AA_simulacoes_lbfgs[,2,,])  # Erro theta
AA_dignostico_lbfgs[,9,,]  <- sqrt(AA_simulacoes_lbfgs[,3,,])  # Erro rho

AA_dignostico_lbfgs[,10,,] <- (( mean(AA_simulacoes_lbfgs[,c(1,2),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_lbfgs[,8,,]) < AA_simulacoes_lbfgs[,5,,]) & (( mean(AA_simulacoes_lbfgs[,c(1,2),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_lbfgs[,8,,]) > AA_simulacoes_lbfgs[,5,,]) # Prob cobertura theta

AA_dignostico_lbfgs[,11,,] <- (( mean(AA_simulacoes_lbfgs[,c(3,4),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_lbfgs[,9,,]) < AA_simulacoes_lbfgs[,6,,]) & (( mean(AA_simulacoes_lbfgs[,c(3,4),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_lbfgs[,9,,]) > AA_simulacoes_lbfgs[,6,,]) # Prob cobertura rho


AA_dignostico_lbfgs[,12,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_lbfgs[,8,,] # Tam_t
AA_dignostico_lbfgs[,13,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_lbfgs[,9,,] # Tam_r
AA_dignostico_lbfgs



save(AA_simulacoes_lbfgs,  file = 'AA_simulacoes_lbfgs.Rdata')
save(AA_dignostico_lbfgs, file = 'AA_diagnostico_lbfgs.Rdata')

AA_simulacoes_lbfgs[,7,1,]
simulacoes_bfgs[,7,,]




AA_simulacoes_lbfgs[,,2,3]


# simulacoes_cg -----------------------------------------------------------

library(LambertW)


par_comb <- list(c(0.5, 0.1), c(0.5, 0.5), c(0.5, 0.8), # Ele pede 9 combinações de parâmetros
                 c(  1, 0.1), c(  1, 0.5), c(  1, 0.8),
                 c(  3, 0.1), c(  3, 0.5), c(  3, 0.8))
length(par_comb)
N <- 50000

AA_simulacoes_cg <- array(c(rep(0,6)), dim=c(N,13,9,10)) # COLUNAS, theta_op1, theta_op2, rho_op1, rho_op2, par1, par2, n, var_theta_op1, var_theta_op2, var_rho_op1, var_rho_op2, cov_theta, cov_rho
AA_simulacoes_cg


set.seed(9999)
for (i in 1:N)                                  # Número de simulações
{ 
  for (index_n in 1:10)                         # Tamanho da amostra
  { n <- seq(10, 100, 10)[index_n]
  
  for (index_par in 1:9)                         # Combinação de parâmetros
  { par <- par_comb[[index_par]]
  
  correlacao <- 1
  while(correlacao > 0)
  { amostra <- rlindley_geom_AA(n, par=par)   # Amostra
  amostra1 <- amostra$amostra1
  amostra2 <- amostra$amostra2
  correlacao <- cor(amostra1, amostra2)}
  
  
  
  op1      <- try(optim(par = par,            # Chute inicial
                        fn = log_lindley_geometrica,  # Log-Verossimilhança
                        x = amostra1,                  # Amostra
                        control = list(fnscale = -1),
                        method = 'CG',              # Método
                        hessian = F                 # Calcular a hessiana
  ))                 
  
  op2      <- try(optim(par = par,            # Chute inicial
                        fn = log_lindley_geometrica,  # Log-Verossimilhança
                        x = amostra2,                  # Amostra
                        control = list(fnscale = -1),
                        method = 'CG',          # Método
                        hessian = F                # Calcular a hessiana
                        
  ))                 
  
  if(typeof(op2) == 'character' || typeof(op1) == 'character')
  { valores <- c(NA, NA, NA, NA, par[1], par[2], n, rep(NA,6))
  AA_simulacoes_cg[i, ,index_par, index_n] <- valores
  next}
  
  valores <- c(op1$par[1], op2$par[1],op1$par[2], op2$par[2], par[1], par[2], n, rep(NA,6))
  # Valores recebe o que queremos dessa bagaça toda,
  # theta_estimado, rho_estimado, theta_real, rho_real, n, se convergiu(0 = sim), variância_rho, variância_theta
  
  cat('itr:', i, '-' , valores, '\n')
  AA_simulacoes_cg[i, ,index_par, index_n] <- valores
  
  }
  }
  
}
AA_simulacoes_cg
load('AA_simulacoes_cg.Rdata')

for(i in 1:10)
{
  for(j in 1:9)
  {
    matriz_theta <- cov(na.omit(AA_simulacoes_cg[,c(1,2),j,i]))
    matriz_rho   <- cov(na.omit(AA_simulacoes_cg[,c(3,4),j,i]))
    
    
    AA_simulacoes_cg[,8,j,i]  <- matriz_theta[1,1]
    AA_simulacoes_cg[,9,j,i]  <- matriz_theta[2,2]
    AA_simulacoes_cg[,10,j,i] <- matriz_rho[1,1]
    AA_simulacoes_cg[,11,j,i] <- matriz_rho[2,2]
    
    AA_simulacoes_cg[,12,j,i] <- matriz_theta[1,2]
    AA_simulacoes_cg[,13,j,i]<- matriz_rho[1,2]
  }
}

AA_dignostico_cg <- array(c(rep(0,6)), dim=c(N,13,9,10))   

AA_dignostico_cg[,1,,]  <- AA_simulacoes_cg[,7,,] # n
AA_dignostico_cg[,2,,]  <- (AA_simulacoes_cg[,8,,] + AA_simulacoes_cg[,9,,])*0.25 + 0.5*AA_simulacoes_cg[,c(12),,] # var theta
AA_dignostico_cg[,3,,]  <-  (AA_simulacoes_cg[,10,,] + AA_simulacoes_cg[,11,,])*0.25 + 0.5*AA_simulacoes_cg[,c(13),,] # var rho
AA_dignostico_cg[,4,,]  <- (AA_simulacoes_cg[,1,,] + AA_simulacoes_cg[,2,,])*0.5 - AA_simulacoes_cg[,5,,]     # Vício de theta
AA_dignostico_cg[,5,,]  <- (AA_simulacoes_cg[,3,,] + AA_simulacoes_cg[,4,,])*0.5 - AA_simulacoes_cg[,6,,]     # Vício de rho
AA_dignostico_cg[,6,,]  <- ((AA_simulacoes_cg[,1,,] + AA_simulacoes_cg[,2,,])*0.5 - AA_simulacoes_cg[,5,,])^2 # EQM de theta
AA_dignostico_cg[,7,,]  <- ((AA_simulacoes_cg[,3,,] + AA_simulacoes_cg[,4,,])*0.5 - AA_simulacoes_cg[,6,,])^2 # EQM de rho
AA_dignostico_cg[,8,,]  <- sqrt(AA_simulacoes_cg[,2,,])  # Erro theta
AA_dignostico_cg[,9,,]  <- sqrt(AA_simulacoes_cg[,3,,])  # Erro rho

AA_dignostico_cg[,10,,] <- (( mean(AA_simulacoes_cg[,c(1,2),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_cg[,8,,]) < AA_simulacoes_cg[,5,,]) & (( mean(AA_simulacoes_cg[,c(1,2),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_cg[,8,,]) > AA_simulacoes_cg[,5,,]) # Prob cobertura theta

AA_dignostico_cg[,11,,] <- (( mean(AA_simulacoes_cg[,c(3,4),,], na.rm=T) - qnorm(1 - 0.05/2) * AA_dignostico_cg[,9,,]) < AA_simulacoes_cg[,6,,]) & (( mean(AA_simulacoes_cg[,c(3,4),,], na.rm=T) + qnorm(1 - 0.05/2) *AA_dignostico_cg[,9,,]) > AA_simulacoes_cg[,6,,]) # Prob cobertura rho


AA_dignostico_cg[,12,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_cg[,8,,] # Tam_t
AA_dignostico_cg[,13,,] <- 2 * qnorm(1 - 0.05/2) * AA_dignostico_cg[,9,,] # Tam_r
AA_dignostico_cg



save(AA_simulacoes_cg,  file = 'AA_simulacoes_cg.Rdata')
save(AA_dignostico_cg, file = 'AA_diagnostico_cg.Rdata')


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



# Resultados --------------------------------------------------------------

# Newton ------------
load('Rdatas/simulacoes_newton.Rdata')
AA_interno_newton <- gera_resultados(AA_dignostico_newton, 'interna', nome_simulacao = 'Com redução', red_var = T)
interno_newton    <- gera_resultados(simulacoes_newton, 'interna', nome_simulacao = 'Sem redução')
resultados_newton_AA <- agrupa_resultados(list(interno_newton, AA_interno_newton))
save(resultados_newton_AA, file ='Rdatas/red_var/resultados/resultados_newton_AA.Rdata')

# Nelder ------------
load('Rdatas/simulacoes_nelder.Rdata')
load('Rdatas/red_var/AA_dignostico_nelder.Rdata')
AA_interno_nelder <- gera_resultados(AA_dignostico_nelder, 'interna', nome_simulacao = 'Com redução', red_var = T)
interno_nelder    <- gera_resultados(simulacoes_nelder, 'interna', nome_simulacao = 'Sem redução')
resultados_nelder_AA <- agrupa_resultados(list(interno_nelder, AA_interno_nelder))
save(resultados_nelder_AA, file ='Rdatas/red_var/resultados/resultados_nelder_AA.Rdata')

# BFGS ------------
load('Rdatas/simulacoes_bfgs.Rdata')
load('Rdatas/red_var/AA_dignostico_bfgs.Rdata')
AA_interno_bfgs <- gera_resultados(AA_dignostico_bfgs, 'interna', nome_simulacao = 'Com redução', red_var = T)
interno_bfgs    <- gera_resultados(simulacoes_bfgs, 'interna', nome_simulacao = 'Sem redução')
resultados_bfgs_AA <- agrupa_resultados(list(interno_bfgs, AA_interno_bfgs))
save(resultados_bfgs_AA, file ='Rdatas/red_var/resultados/resultados_bfgs_AA.Rdata')

# LBFGS ------------
load('Rdatas/simulacoes_lbfgs.Rdata')
AA_interno_lbfgs <- gera_resultados(AA_dignostico_lbfgs, 'interna', nome_simulacao = 'Com redução', red_var = T)
interno_lbfgs    <- gera_resultados(simulacoes_lbfgs, 'interna', nome_simulacao = 'Sem redução')
resultados_lbfgs_AA <- agrupa_resultados(list(interno_lbfgs, AA_interno_lbfgs))
save(resultados_lbfgs_AA, file ='Rdatas/red_var/resultados/resultados_lbfgs_AA.Rdata')

# CG ---------------
load('Rdatas/simulacoes_cg.Rdata')
AA_interno_cg <- gera_resultados(AA_dignostico_cg, 'interna', nome_simulacao = 'Com redução', red_var = T)
interno_cg    <- gera_resultados(simulacoes_cg, 'interna', nome_simulacao = 'Sem redução')
resultados_cg_AA <- agrupa_resultados(list(interno_cg, AA_interno_cg))
save(resultados_cg_AA, file ='Rdatas/red_var/resultados/resultados_cg_AA.Rdata')

# SANN -------------
load('Rdatas/simulacoes_sann.Rdata')
load('Rdatas/red_var/AA_diagnostico_sann.Rdata')
AA_interno_sann <- gera_resultados(AA_diagnostico_sann.Rdata, 'interna', nome_simulacao = 'Com redução', red_var = T)
interno_sann    <- gera_resultados(simulacoes_sann.Rdata, 'interna', nome_simulacao = 'Sem redução')
resultados_sann_AA <- agrupa_resultados(list(interno_sann, AA_interno_sann))
save(resultados_sann_AA, file ='Rdatas/red_var/resultados/resultados_sann_AA.Rdata')



# gráficos ----------------------------------------------------------------

library(ggplot2)

formato <- theme(                                                       
  plot.title = element_text(size = 14, hjust = 0.5),
  axis.title.y = element_text(size = 12, vjust = 0.5, angle= 0),
  axis.title.x = element_text(size = 12, vjust = -0.2),
  axis.text.y = element_text(size = 10),
  axis.text.x = element_text(size = 10)
)
resultados_bfgs_AA



reducao_lineplot <- function(var, banco, simulacao, xlab = 'Tamanho da amostra', paleta = NULL)
{
  
  
  
  ifelse(var == 'var_t'   || var == 'var_r',  titulo <- 'Variância média para',
         ifelse(var == 'vicio_t' || var == 'vicio_r',titulo <- 'Vício médio para',
                ifelse(var == 'eqm_t'   || var == 'eqm_r',  titulo <- 'EQM médio para',
                       ifelse(var == 'erro_t'  || var == 'erro_r', titulo <- 'Erro padrão médio para',
                              ifelse(var == 'prob_t'  || var == 'prob_r', titulo <- 'Probabilidade média da cobertura para',
                                     ifelse(var == 'tam_t'   || var == 'tam_r',  titulo <- 'Tamanho médio da cobertura para'))))))
  
  ylab <- sub(' .*', '',titulo)
  ifelse(substr(var,nchar(var),nchar(var)) == 'r',
         titulo <- as.expression(bquote(.(titulo) ~ rho)),
         titulo <- as.expression(bquote(.(titulo) ~ theta)))  
  
  
  reducao <- unique(banco$simulacao)
  banco1  <- banco[banco$simulacao == reducao[1],]
  banco2  <- banco[banco$simulacao == reducao[2],]
  
  if(is.null(paleta)) paleta <- c('#004586', '#FF420E', '#FFD320', '#579D1C', '#7E0021', '#83CAFF', '#FF950E', '#AECF00', '#DD4477')     
  
  if(var == 'prob_t' || var == 'prob_r') intercept <- 0.95
  else intercept <- 0
  
  grafico1 <- ggplot(data = banco1, aes(x = n, y = as.numeric(.data[[var]]), colour = comb_par, group = comb_par)) +
    geom_line( linewidth =0.9, alpha = 0.7, linetype = 'dashed')+
    geom_point(size=2)+
    xlab(xlab)+
    ylab(ylab)+
    ggtitle(titulo)+
    theme_minimal()+
    formato +
    scale_colour_manual(values=paleta)+
    guides(colour = guide_legend(title =  ~~~~~~~~~ theta ~~~~ rho))+
    geom_hline(yintercept=intercept, col = 'red', linetype = "dotdash")
  
  grafico2 <- ggplot(data = banco2, aes(x = n, y = as.numeric(.data[[var]]), colour = comb_par, group = comb_par)) +
    geom_line( linewidth =0.9, alpha = 0.7, linetype = 'dashed')+
    geom_point(size=2)+
    xlab(xlab)+
    ylab(ylab)+
    ggtitle(titulo)+
    theme_minimal()+
    formato +
    scale_colour_manual(values=paleta)+
    guides(colour = guide_legend(title =  ~~~~~~~~~ theta ~~~~ rho))+
    geom_hline(yintercept=intercept, col = 'red', linetype = "dotdash")
  
  return(grafico)
}



graficos_reducao <- function(banco, simulacao)
{ 
  variaveis_diag <- c('var_t', 'var_r', 'eqm_t', 'eqm_r','vicio_t','vicio_r', 'erro_t', 'erro_r', 'prob_t', 'prob_r','tam_t', 'tam_r')

    for(i in 1:(length(variaveis_diag)/3))
    {
      plots1 <- comb_par_lineplot(variaveis_diag[i*4 -3], banco, simulacao)
      plots2 <- comb_par_lineplot(variaveis_diag[i*4 -2], banco, simulacao)
      plots3 <- comb_par_lineplot(variaveis_diag[i*4 -1], banco, simulacao)
      plots4 <- comb_par_lineplot(variaveis_diag[i*4   ], banco, simulacao)
      
      plot_final    <- ggarrange(plots1, plots2,plots3, plots4, nrow = 4, common.legend = T, legend = 'right')
      
      arquivo <- paste0('img/red_var/',simulacao, '/', i , '.jpg')
      ggsave(arquivo, plot = plot_final)
   }
  
}
library('ggplot2')

graficos_reducao(resultados_bfgs_AA, 'BFGS')
