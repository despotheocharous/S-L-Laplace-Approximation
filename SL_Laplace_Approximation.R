#####################################
#
# S-L Laplace Approximation
# 
# Spiegelhalter and Lauritzen (1990)
#
#####################################


SL_Laplace_Approximation <- function(mu_prior, Sigma_prior, X, Y, tol){
  
  # @param mu_prior: prior mean vector 
  # @param Sigma_prior: prior Sigma matrix
  # @param X: design matrix n x (p+1) (including the intercept column)
  # @param Y: observation vector
  # @param tol: tolerance for convergence
  # @return mu_post: posterior mean vector
  # @return Sigma_post: posterior Sigma matrix
  # @return iterations: number of iterations needed for convergence
  
  iter <- 0
  
  mu_old <- mu_prior
  Sigma_old <- Sigma_prior
  
  continue_indicator <- TRUE
  
  while (continue_indicator){
    
    iter <- iter + 1
    
    theta <- X %*% mu_old
    pi <- 1/(1+exp(-theta))
    pi_ <- pi * (1-pi) # pi element-wise multiplication with 1-pi 
    
    Sigma_post <- solve(solve(Sigma_old) + t(X)%*%diag(as.list(pi_))%*%X)
    mu_post <- mu_old + Sigma_post%*%t(X)%*%(Y-pi)
    
    # mu and Sigma must converge
    convergence_condition1 <- all(abs(Sigma_post-Sigma_old) < tol)
    convergence_condition2 <- all(abs(mu_post-mu_old) < tol)
    
    if (convergence_condition1 && convergence_condition2) {
      continue_indicator <- FALSE
    }
    
    # Udpate mu and Sigma for the next iteration
    mu_old <- mu_post
    Sigma_old <- Sigma_post
  }
  
  return(list(mu_post = mu_post, 
              Sigma_post = Sigma_post, 
              iterations = iter))
}