normalise_logprob <- function(x){
  x <- exp(x-max(x))
  return(x/sum(x))
}

logsumexp <- function(x){
  log(sum(exp(x - max(x)))) + max(x)
}

# Corrected: removed n,p
get_clusterprobs <- function(dat, g, distr, mu, sigma, dof=NULL, delta=NULL, pro){
  logdens <- ddmix(dat=dat, g=g, distr=distr, mu=mu, sigma=sigma, dof=dof, delta=delta)
  logprobs <- t(t(logdens)+log(pro))
  clusprobs <- t(apply(logprobs, 1, normalise_logprob))
  return(clusprobs)
}

# Corrected: removed n,p
get_entropy <- function(dat, g, distr, mu, sigma, dof=NULL, delta=NULL, pro){
  tau <- get_clusterprobs(dat, g, distr, mu, sigma, dof=NULL, delta=NULL, pro)
  entropy <- apply(tau, 1, function(x) sum(ifelse(x==0, 0,-x*log(x))))
  return(entropy)
}

makelabelmatrix <- function(clust){
  n <- length(clust)
  g <- max(clust)
  Z <- matrix(0, n, g)
  index <- cbind(1:n, clust)
  Z[index] <- 1
  return(Z)
}


lrlk <- function(y, X, beta){
  eta <- X %*% beta
  loglk <- sum(beta*(t(X) %*% y))-sum(sapply(eta, function(z) logsumexp(c(0, z))))
  return(loglk)
}


cov2vec <- function(sigma){
  R <- chol(sigma)
  upper_elements <- R[upper.tri(R, diag=FALSE)]
  diag_elements <- diag(R)
  return(c(log(diag_elements), upper_elements))
}


vec2cov <- function(par){
  q <- length(par)
  p <- (-1+sqrt(1+4*1*2*q))/2
  R <- matrix(0, p,p)
  diag_elements <- par[1:p]
  upper_elements <- par[-(1:p)]
  diag(R) <- exp(diag_elements)
  if (any(is.infinite(R))){
    stop('Variances infinite')
  }
  if (any(is.nan(R))){
    #browser() 
    stop('Variances NaN')
  }
  R[upper.tri(R)] <- upper_elements
  sigma <- t(R) %*% R
  if(any(eigen(sigma)$values <= 0)) stop('Conversion failed (negative eigenvalues)')
  return(sigma)
}

pro2vec <- function(pro){
  g <- length(pro)
  z <- numeric(g)
  z[1] <- pro[1]
  for (h in 2:g){
    z[h] <- pro[h]/(1-sum(pro[1:(h-1)]))
  }
  z <- z[-g]
  y <- log(z/(1-z))-log(1/(g-1:(g-1)))
  return(y)
}

vec2pro <- function(vec){
  g <- length(vec)+1
  pro <- numeric(g)
  z <- numeric(g)
  for (h in 1:(g-1)){
    z[h] <- 1/(1+exp(-vec[h]-log(1/(g-h))))
  }
  pro[1] <- z[1]
  if (g > 2){
    for (h in 2:(g-1)){
      pro[h] <- (1-sum(pro[1:(h-1)]))*z[h]
    }
  }
  pro[g] <- 1-sum(pro[-g])
  return(pro)
}

par2list <- function(par, g, p){
  mu <- matrix(par[1:(p*g)], p, g)
  par <- par[-(1:(p*g))]
  q <- p*(p+1)/2
  cholpars <- matrix(par[1:(q*g)], q, g)
  sigma <- array(0, dim=c(p,p,g))
  for (h in 1:g){
    sigma[,,h] <- vec2cov(cholpars[,h])
  }
  par <- par[-(1:(q*g))]
  tpro <- par[1:(g-1)]
  par <- par[-(1:(g-1))]
  parlist <- list(mu=mu, sigma=sigma, pro=vec2pro(tpro))
  return(parlist)
}

# Cholesky decomposition and vectorisation of sigma (via cov2vec),
# then 
list2par <- function(parlist, g, p){
  muvec <- as.vector(parlist$mu)
  q <- p*(p+1)/2
  cholpars <- matrix(0, q, g)
  for (h in 1:g){
    cholpars[,h] <- cov2vec(parlist$sigma[,,h]) # Cholesky decomposition and vec
  }
  cholvec <- as.vector(cholpars)
  tpro <- pro2vec(parlist$pro)
  par <- c(muvec, cholvec, tpro)
  return(par)
}


neg_profile_objective_function <- function(par, X1, Z1, X2, m, g){
  # for use with nlminb (function for minimisation)
  -profile_objective_function(par, X1, Z1, X2, m, g)
}

profile_objective_function <- function(par, X1, Z1, X2, m, g){
  # returns profile likelihood, given estimates of the mixture parameters,
  # ml estimates of xi_0 and xi_1 are given by glm function
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  p <- ncol(X1)
  
  parlist <- par2list(par, g, p)
  
  mu <- parlist$mu
  sigma <- parlist$sigma
  pro <- parlist$pro

  
  dat <- rbind((X1), (X2))
  
  tau <- get_clusterprobs(dat=X2, g=g, distr='mvn', mu=mu, sigma=sigma, delta=NULL, pro=pro)
  unlabelled <- ddmix(dat=X2, g=g, distr='mvn', mu=mu, sigma=sigma, delta=NULL)
  labelled <- ddmix(dat=X1, g=g, distr='mvn', mu=mu, sigma=sigma, delta=NULL)
  
  # Missingness mechanism log likelihood assuming eta is simple linear function of entropy
  entropy <- get_entropy(dat, g, distr='mvn', mu=mu, sigma=sigma, pro=pro)
  # design matrix for missingness model
  D <- cbind(1, -(entropy))
  glmfit <- glm(m ~ 0 + D, family=binomial)
  beta <- coef(glmfit)
  mislk <- lrlk(m, D, beta)
  
  
  unlablk <- sum(log(colSums(exp(t(unlabelled))*pro)))
  lablk <- sum(log(colSums(exp(t(labelled))*pro*t(Z1))))

  value <- mislk + unlablk + lablk
  return(value)
}


