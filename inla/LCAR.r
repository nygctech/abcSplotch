'inla.rgeneric.LCAR.model' <- function(cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
                                               "log.prior", "quit"), theta = NULL){
  interpret.theta <- function()
  {
    alpha <- 1 / (1 + exp(-theta[1L]))
    prec <- exp(theta[2L])
    param = c(alpha, prec)
    return(list(alpha = alpha, prec = prec, param = param))
  }
  
  graph <- function()
  {
    G <- Matrix::Diagonal(nrow(W), 1) + W
    return (G)
  }
  
  #Precision matrix
  Q <- function()
  {
    #Parameters in model scale
    param <- interpret.theta()
    n <- nrow(W)

    # Sparse degree vector (# neighbors for each node) -- avoid densification
    d <- Matrix::rowSums(W)

    I <- Matrix::Diagonal(n, x = 1)
    D <- Matrix::Diagonal(n, x = d)
    
    #Precision matrix
    Q <- param$prec * (param$alpha * I + (1 - param$alpha) * (D - W))
    return (Q)
  }
  
  mu <- function() {
    return(numeric(0))
  }
  
  log.norm.const <- function() {
    ## return the log(normalising constant) for the model
    #param = interpret.theta()
    #
    #val = n * (- 0.5 * log(2*pi) + 0.5 * log(prec.innovation)) +
    #0.5 * log(1.0 - param$alpha^2)
    
    val <- numeric(0)
    return (val)
  }
  
  log.prior <- function() {
    ## return the log-prior for the hyperparameters
    param <- interpret.theta()
    
    # log-Prior for the autocorrelation parameter
    val <- - theta[1L] - 2 * log(1 + exp(-theta[1L]))
    
    # Chisquare for precision
    val <- val + dchisq(param$prec, 1, log=T) + log(param$prec)
    return (val)
  }
  
  initial <- function() {
    ## return initial values
    return (c(0,0))
  }
  
  quit <- function() {
    return (invisible())
  }
  
  val <- do.call(match.arg(cmd), args = list())
  return (val)
}

inla.LCAR.model <- function(...) {
  INLA::inla.rgeneric.define(inla.rgeneric.LCAR.model, ...)
}