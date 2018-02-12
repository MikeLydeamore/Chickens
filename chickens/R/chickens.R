#Test chickens model

betas <- matrix(1.5, dimnames=list(c("Sc")))
#Scavenging:
x0 <- list("Sc.E"=1000, "Sc.Ch.S"=100)
n <- c(10, 10, 10)
delta <- c(10, 10, 10, 10)
rho <- 0.5
y <- 0.3
x <- 0.3
alpha <- list("lG"=4)
sigma <- 1
gamma <- 3

p_sub <- list("x0"=x0, "n"=n, "delta"=delta, "rho"=rho, "y"=y, "x"=x, "alpha"=alpha, "sigma"=sigma, "gamma"=gamma)
p <- list("Sc"=p_sub)


runChickensModel <- function(parameter_list, betas)
{
  num_patches <- length(parameter_list)
  
  for (i in 1:num_patches)
  {
    params <- parameter_list[[i]]
    if (!exists("x0", where = params))
      stop("Patch labelled ",names(parameter_list)[i], " requires x0.")
    if (!exists("n", where = params))
      #TODO
      params$n <- c(10, 10, 10)
    if (!exists("delta", where = params))
      params$delta <- c(10, 10, 10, 10)
    if (!exists("y"), where = params)
      params$y <- 0
    if (y < 0 || y > 1)
      stop("y must be between 0 and 1 (inclusive)")
    if (!exists("x"), where = params)
      params$x <- 0
    if (x < 0 || x > 1)
      stop("x must be between 0 and 1 (inclusive")
    if (!exists("alpha", where = params))
      params$alpha <- list("none"=0)
    if (!exists("sigma", where = params))
      params$sigma <- 1
    if (!exists("gamma", where = params))
      params$gamma <- 3
    
    #Gotta feed it back in because R won't do references :'-(
    parameter_list[[i]] <- params
  }
  
  run <- as.data.frame(.chickens_model(parameter_list, betas))
  
  return (run)
}

