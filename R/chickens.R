ggplotMatrixLinesOnly <- function (df, id.vars, ...) 
{
  foo <- melt(df, id.vars)
  geom_line(data = foo, aes_string(x = id.vars, y = "value", 
                                   colour = "variable"), ...)
}

#' Run Chicken Model
#' 
#' Runs a single realisation of the chickens model
#' 
#' @param parameter_list A list of parameters for this realisation. Needs the following structure:
#'   \itemize{
#'   \item{\code{patchname}: One of "Es","Ns","Bs" or "Sc"}
#'     \itemize{
#'        \item{\code{x0}: Initial condition}
#'        \item{\code{delta}: Numeric vector of length 5}
#'        \item{\code{y}: Numeric in [0, 1]}
#'        \item{\code{x}: Numeric in [0, 1]}
#'        \item{\code{alpha}: List containing 3 elements, "lG","He" and "Rs"}
#'        \item{\code{sigma}:}
#'        \item{\code{gamma}:}
#'        \item{\code{w}:}
#'        \item{\code{n_egg}:}
#'        \item{\code{K}: Carrying capacity (Sc system only)}
#'        \item{\code{q}:}
#'     }
#'   }
# 
#'   Repeat for each possible patch.
#'
#' @param betas Matrix of within and between patch transmission (row names are required)
#' @param dt Time spacing of outputs (NOT solving points)
#' @param max_time Maximum time for simulation
#' @param solver_type Either "stochastic" or "deterministic"
#' 
#' @examples 
#' betas <- matrix(1.5, dimnames=list(c("Es")))
#' x0 <- list("E"=100, "Ch.S"=50, "He.S"=50, "He.I"=10)
#' p_sub <- list("x0"=x0)
#' p <- list("Es"=p_sub)
#' df <- runChickensModel(parameter_list = p, betas=betas)
#' 
#' @return A list containing two elements: \code{realisation}, contains the realisation and \code{parameters} contains the parameters
runChickensModel <- function(parameter_list, betas = matrix(), dt = 1, max_time = 1000, solver_type = "stochastic", seed = -1)
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
      params$alpha <- list("lG"=0.1, "He"=0.1, "Rs"=0.1)
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
        params$n_egg <- 8.9
      if (patch_type == "Ns")
        params$n_egg <- 10.8
      if (patch_type == "Bs")
        params$n_egg <- 11.3
      if (patch_type == "Es")
        params$n_egg <- 10.2
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

    if (!exists("br", where=params))
    {
      if (patch_type == "Ns")
        params$br <- 1
      if (patch_type == "Bs")
        params$br <- 1
      if (patch_type == "Es")
        params$br <- 0.1
    }
    
    #Gotta feed it back in because R won't do references :'-(
    parameter_list[[i]] <- params
  }
  if (solver_type == "stochastic")
    solver <- -1
  else
    solver <- 1
  run <- as.data.frame(.chickens_model(parameter_list, betas, max_time, dt, solver, seed))
  
  return (list("realisation"=run, "parameters"=parameter_list, "seed"=seed))
}

#' Get number of chickens at given time
#' 
#' @param state State vector
#' @return Number of chickens (total) in the given state
getNumberOfChickens <- function(state)
{
  cols <- colnames(state)
  idx <- sapply(strsplit(cols, "\\."), length) == 3
  return (sum(state[idx]))
}

#' Get the number of chickens at all times
#' 
#' @param realisation Realisation from the Chickens model
#' @return Vector consisting of the number of chickens at each time
#' 
#' @examples 
#' df <- runChickensModel(parameter_list = p, betas=betas)
#' df$realisation$total_chickens <- getNumberOfChickensAtAllTimes(df$realisation)
getNumberOfChickensAtAllTimes <- function(realisation)
{
  sapply(1:nrow(realisation), function(i) {getNumberOfChickens(realisation[i,])})
}