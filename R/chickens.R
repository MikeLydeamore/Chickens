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


runChickensModel <- function(parameter_list, betas = matrix(), dt = 1, max_time = 1000, solver_type == "stochastic")
{
  num_patches <- length(parameter_list)
  
  for (i in 1:num_patches)
  {
    params <- parameter_list[[i]]
    patch_type <- names(parameter_list)[i]
    if (!exists("x0", where = params))
      stop("Patch labelled ",names(parameter_list)[i], " requires x0.")
    if (!exists("n", where = params))
      params$n <- c(1/21, 1/28, 1/70)
    if (!exists("delta", where = params))
    {
      if (patch_type == "Sc")
        params$delta <- c(1/47.05, 1/25.26, 1/167.79, 1/13718, 1/4801.3)
      if (patch_type == "Ns")
        params$delta <- c(1/28.05, 1/107.1, 1/332.2,  1/480.1, 1/480.1)
      if (patch_type == "Bs")
        params$delta <- c(1/28.05, 1/125.5, 1/322.2, 1/480.1, 1/480.1)
      if (patch_type == "Es")
        params$delta <- c(1/28.05, 1/172.3, 1/502.6, 1/480.1, 1/480.1)
    }
    if (!exists("y", where = params))
    {
      if (patch_type == "Sc")
        params$y <- 0
      if (patch_type == "Ns")
        params$y <- 1
      if (patch_type == "Bs")
        params$y <- 1
      if (patch_type == "Es")
        params$y <- 0.8
    }
    if (params$y < 0 || params$y > 1)
      stop("y must be between 0 and 1 (inclusive)")
    if (!exists("x", where = params))
    {
      if (patch_type == "Sc")
        params$x <- 0.83
      if (patch_type == "Ns")
        params$x <- 0.60
      if (patch_type == "Bs")
        params$x <- 0.75
      if (patch_type == "Es")
        params$x <- 0.62
    }
    if (params$x < 0 || params$x > 1)
      stop("x must be between 0 and 1 (inclusive)")
    if (!exists("alpha", where = params))
      params$alpha <- list("lG"=1, "He"=1, "Rs"=3)
    if (!exists("sigma", where = params))
      params$sigma <- 1
    if (!exists("gamma", where = params))
      params$gamma <- 3
    
    if (!exists("w", where = params))
    {
      if (patch_type == "Sc")
        params$w <- 1/32.1
      if (patch_type == "Ns")
        params$w <- 1/2.6
      if (patch_type == "Bs")
        params$w <- 1/2.7
      if (patch_type == "Es")
        params$w <- 1/2.6
    }
    
    if (!exists("n_egg", where=params))
    {
      if (patch_type == "Sc")
        params$n_egg <- 1/32.1
      if (patch_type == "Ns")
        params$n_egg <- 1/2.6
      if (patch_type == "Bs")
        params$n_egg <- 1/2.7
      if (patch_type == "Es")
        params$n_egg <- 1/2.6
    }
    
    if (!exists("K", where=params))
      params$K <- 1100
    
    if (!exists("q", where=params))
    {
      if (patch_type == "Sc")
        params$q <- 4
      if (patch_type == "Ns")
        params$q <- 2
      if (patch_type == "Bs")
        params$q <- 2
      if (patch_type == "Es")
        params$q <- 2
    }
    #Gotta feed it back in because R won't do references :'-(
    parameter_list[[i]] <- params
  }
  if (solver_type == "stochastic")
    solver <- -1
  else
    solver <- 1
  run <- as.data.frame(.chickens_model(parameter_list, betas, dt, solver))
  
  return (run)
}