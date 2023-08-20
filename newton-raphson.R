dtheta1 <- (((2*x^2+2*x)*(theta^2)+(2*x^2+4*x)*theta)*rho)
dtheta2 <- (((x+1)*theta^2+(x+2)*theta+1)*rho+(-(theta^2)-2*theta-1)*exp(x*theta))
dtheta  <- 2*n/theta - n/(theta + 1) - sum(x) + sum(dtheta1/dtheta2)


par[1] = 2

par[2] = 0.5

x = rlingley_geom(1000, c(2, 0.5))
n = length(x)
-sum(((2*(x*((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2]-(x/(par[1]+1)-(x*par[1])/(par[1]+1)^2)*exp(-(x*par[1]))*par[2]))/(1-((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2])))-n/(par[1]+1)+(2*n)/par[1]-sum(x)



U         <- function(x, par)
{    # Vetor escore é um vetor que contém a primeira derivada em relação a cada parâmetro, nesse caso é a derivada primeira em relação a theta, derivada primeira em relação a rho
  theta <- par[1]
  rho   <- par[2]
  n       <- length(x)
  
  dtheta1 <- (((2*x^2+2*x)*(par[1]^2)+(2*x^2+4*x)*par[1])*par[2])
  dtheta2 <- (((x+1)*par[1]^2+(x+2)*par[1]+1)*par[2]+(-(par[1]^2)-2*par[1]-1)*exp(x*par[1]))
  dtheta  <- 		sum((2*(x*((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2]-(x/(par[1]+1)-(x*par[1])/(par[1]+1)^2)*exp(-(x*par[1]))*par[2])^2)/(1-((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2])^2)
    -sum((2*(-(x^2*((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2]))))
        +sum(2*x*(x/(par[1]+1)-(x*par[1])/(par[1]+1)^2)*exp(-(x*par[1]))*par[2])
        -sum(((2*x*par[1])/(par[1]+1)^3-(2*x)/(par[1]+1)^2)*exp(-(x*par[1]))*par[2]/(1-((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2]))
  +n/(par[1]+1)^2-(2*n)/par[1]^2
  drho1   <- ((x*theta)/(theta+1)+1)
  drho2   <- exp(-(x*theta))
  drho    <- sum((2*((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1])))/(1-((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2]))-n/(1-par[2])
  
  vetor   <- matrix(nrow = 2, c(dtheta, drho))
  
  return(vetor)
}

set.seed(123)
U(x, c(2,0.3))




dtheta2   <- function(x, par)           # Derivada segunda da log-verossimilhança em relação a theta
{
  n       <- length(x)
  f1      <- (x*par[1])/(par[1] + 1) + 1
  f2      <- exp(-x*par[1])*par[2]
  f3      <- x/(par[1] + 1)
  f4      <- (x*par[1])/((par[1] + 1)^2)
  
  d1      <- 2*(((x*(f1)*f2) - (f3 - f4)*f2)^2)/(1 - (f1)*f2)^2
  d2      <- (2*(-(x^2*(f1)*f2) + 2*x*(f3 - f4)*f2 - ((2*x*par[1])/(par[1] + 1)^3 - 2*x/(par[1] + 1)^2)*f2))/(1 - f1*f2)
  
  dd      <- sum((2*(x*((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2]-(x/(par[1]+1)-(x*par[1])/(par[1]+1)^2)*exp(-(x*par[1]))*par[2])^2)/(1-((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2])^2)-sum(
    (2*(-(x^2*((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2]))+2*x*(x/(par[1]+1)-(x*par[1])/(par[1]+1)^2)*exp(-(x*par[1]))*par[2]))-
    sum(((2*x*par[1])/(par[1]+1)^3-(2*x)/(par[1]+1)^2*exp(-(x*par[1]))*par[2])/(1-((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2]))+n/(par[1]+1)^2-(2*n)/par[1]^2
  
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
  
  dd      <- sum((2*((x*par[1])/(par[1]+1)+1)^2*exp(-(2*x*par[1])))/(1-((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2])^2)-n/(1-par[2])^2
  
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
  
  dd <-  -sum(((2*(x*((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))-(x/(par[1]+1)-(x*par[1])/(par[1]+1)^2)*exp(-(x*par[1]))))/(1-((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2])))-sum((2*((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*(x*((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2]-(x/(par[1]+1)-(x*par[1])/(par[1]+1)^2)*exp(-(x*par[1]))*par[2]))/(1-((x*par[1])/(par[1]+1)+1)*exp(-(x*par[1]))*par[2])^2)
  
  
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

J(rlingley_geom(1000, c(0.3, 0.5)), c(0.3, 0.5))





par <- matrix(c(2, 0.1), nrow = 2)                # Os valores que o Newton-raphson precisa
set.seed(9999)
x <-  rlingley_geom(50000, c(2, 0.1))
itr <- 5000
erro <- 0.0001
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

a <- c(rep(0, 5000))
for (i in i:5000)
{
 a[i] <- log_lindley_geometrica(rlingley_geom(5000, c(2, 0.1)), c(2,0.1))
}
plot(a)
