
treeAgr.EM <-
  function(S,
           k,
           K0 = NULL,
           Sigma0 = NULL,
           pii=0.5,
           n,
           max.iter = 100,
           eps = 0.1) {
    # -----------------------------------------------------------------------------------------------------------------------------
    # FUNCTION
    #   Perform EM with Matrix Tree theorem when k variable are missing
    # INPUT
    #   S : empirical covariance matrix of observed variables (pxp matrix)
    #   Sigma0 : initial value of complete covariance matrix ((p+k)x(p+k) matrix)
    #   k : number of hidden variables (integer)
    #   max.iter : number of iterations (integer, default is 10)
    # OUTPUT
    #   Sigma    : estimator of complete covariance matrix ((p+k)x(p+k) matrix)
    #   Alpha    : posterior edge probabilities ((p+k)x(p+k) matrix)
    #   K        : expectation of precision matrix under posterior edge distribution ((p+k)x(p+k) matrix)
    #   ordered_couples : list of edges in the decreasing order of the posterior probabilities
    #     list of size p(p-1)/2. The k_th element  of the list noted e is a vector of size 2 such that
    #     i = e[1] and j=e[2] are the indices of the nodes that have the k_th highest posterior probability
    #     of being connected by an edge
    #   likelihood: vector containing marginal likelihood at every iteration
    # -----------------------------------------------------------------------------------------------------------------------------

    p <- ncol(S)
    r <- ncol(Sigma0) - p
    if (r == 0) {
      Sigma=S
      K = solve(S)
      alpha = matrix(rep(pii, p*p), p)
      likelihoods = rep(0,max.iter)
      likelihoods[1] = -log(det(K))+ matrix.trace(K %*% S)-n*log(2*pi)/2
      error = rep(0,max.iter)
      error_ = rep(0,max.iter)
      error[1]=1
      error_[1]=1
      it=1
      while (error_[it] > eps && it<=max.iter) {
        print(paste('EM iteration', it, sep = " "))
        alpha.tmp=alpha
        D <- diag(K) %*% t(diag(K))
        d <- ((D - K ^ 2) / D)
        diag(d) = 1
        d = log(d)
        diag(d) = log(diag(K))
        t <- -2 * K * S
        diag(t) = -diag(K) * diag(S)
        gamma <- d + t
        diag(gamma) = 0
        magic<-1/n
        alpha <- abs(edge.prob(magic*n*gamma,account.prior = TRUE, q0 = pii))
        #while (sum(is.na(alpha))>0){
        #  magic<-magic/2
        #  alpha <- abs(edge.prob(magic*n*gamma,account.prior = TRUE, q0 = pii))
        #}
        alpha[alpha>1]=1
        alpha[is.na(alpha)]=0
        ## M step:
        K2 <- matrix(0, p, p)
        wrap_optimDiag <- function(k_ii, i) {
          return(optimDiag(i, k_ii, K, S, alpha))
        }
        maxi = 1e50
        mini = 1e-50
        K2diagonal <-
          sapply(1:p, function(i)
            dichotomie(mini, maxi, function(x)
              wrap_optimDiag(x, i), 1e-5))
        diag(K2) = K2diagonal
        D2 = diag(K2) %*% t(diag(K2))
        K2 <- (1 - sqrt(1 + 4 * S ^ 2 * D2[1:p, 1:p])) / (2 * S)
        diag(K2) = diag(K)
        K1 <- nearPD(K2, eig.tol = 0.1)
        K <- as.matrix(K1$mat)
        likelihoods[it] <- -log(det(K)) + matrix.trace(K %*% S)+log(2*pi)/2
        error[it] = norm(alpha - alpha.tmp,type="F")
        it=it+1
        error_[it] = abs(likelihoods[it] - likelihoods[it - 1])
      }
    }else{
      Sigma <- Sigma0
      K <- K0
      it <- 1
      if (k == 1) {
        Km <-
          K[1:p, 1:p] - (1 / K[p + 1, p + 1]) * K[1:p, (p + 1):(p + r)] %*% t(K[(p +
                                                                                   1):(p + r), 1:p])
      } else {
        Km <-
          K[1:p, 1:p] - K[1:p, (p + 1):(p + r)] %*% solve(K[(p + 1):(p + r), (p +
                                                                                1):(p + r)]) %*% K[(p + 1):(p + r), 1:p]
      }
      likelihoods <- rep(0, max.iter + 1)
      error <- rep(0, max.iter + 1)
      error_ <- rep(0, max.iter + 1)
      alpha = matrix(rep(pii, (p + r) * (p + r)), p + r)
      gamma = log(matrix(rep(pii, (p + r) * (p + r)), p + r))/n
      P = log(matrix(rep(pii, (p + r) * (p + r)), p + r))
      B = log(matrix(rep(pii, r * r), r))
      W = log(matrix(rep(pii, r * p), r))
      diag(alpha) = 0
      it = 1
      error[it] = 1
      error_[it] = 1
      likelihoods[it] <- -log(det(Km)) + matrix.trace(Km %*% S)-log(2*pii)/2
      while (it <= max.iter && error_[it] > eps) {
        print(paste('EM iteration', it, sep = " "))

        ## E step:
        Sigma.tmp <- Sigma
        Sigma[(p + 1):(p + r), 1:p] <-
          -solve(K[(p + 1):(p + r), (p + 1):(p + r)]) %*% K[(p + 1):(p + r), 1:p] %*% Sigma[1:p, 1:p]
        Sigma[1:p, (p + 1):(p + r)] <-
          -Sigma[1:p, 1:p] %*% K[1:p, (p + 1):(p + r)] %*% solve(K[(p + 1):(p + r), (p + 1):(p + r)])
        Sigma[(p + 1):(p + r), (p + 1):(p + r)] <-
          solve(K[(p + 1):(p + r), (p + 1):(p + r)]) + solve(K[(p + 1):(p + r), (p + 1):(p + r)]) %*%
          K[(p + 1):(p + r), 1:p] %*% Sigma[1:p, 1:p] %*% K[1:p, (p + 1):(p + r)] %*%
          solve(K[(p + 1):(p + r), (p + 1):(p + r)])

        D <- diag(K) %*% t(diag(K))
        d <- ((D - K ^ 2) / D)
        diag(d) = 1
        d = log(d)
        diag(d) = log(diag(K))
        t <- -2 * K[1:p, 1:p] * Sigma[1:p, 1:p]
        diag(t) = -diag(K[1:p, 1:p]) * diag(Sigma[1:p, 1:p])
        f <- computeF(K, Sigma[1:p, 1:p], p)
        m <- matrix(0, p + r, p + r)
        m[1:p, 1:p] <- t
        m[1:p, (p + 1):(p + r)] <- f
        m[(p + 1):(p + r), 1:p] <- t(f)
        gamma <- d + m + log(matrix(rep(pii, (p + r) * (p + r)), p + r))/n
        diag(gamma) = 0
        alpha.tmp = alpha
        magic<-1/n
        alpha <- abs(edge.prob(magic*n*gamma,account.prior = TRUE, q0 = pii))
        #while (sum(is.na(alpha))>0){
        #  magic<-magic/2
        #  alpha <- abs(edge.prob(magic*n*gamma,account.prior = TRUE, q0 = pii))
        #}
        alpha[alpha>1]=1
        alpha[is.na(alpha)]=0
        W <- solve(K[(p + 1):(p + r), (p + 1):(p + r)]) %*%K[(p + 1):(p + r), 1:p] %*% S
        V = solve(K[(p + 1):(p + r), (p + 1):(p + r)]) %*% K[(p + 1):(p + r), 1:p]
        V = V %*% Sigma[1:p, 1:p] %*% t(V)
        B <- solve(K[(p + 1):(p + r), (p + 1):(p + r)]) + V


        ## M step:
        K2 <- matrix(0, p + r, p + r)
        wrap_optimDiag <- function(k_ii, i) {
          return(optimDiag(i, k_ii, K, Sigma, alpha))
        }
        wrap_optimDiagH <- function(k_ii, i) {
          return(optimDiag(i, k_ii, K, B, alpha))
        }

        maxi = 1e50
        mini=1e-50
        #mini = sapply(1:(p+r), function(idx) max(K[idx,setdiff(1:(p+r),idx)]^2/diag(K)[setdiff(1:(p+r),idx)]))
        #mins <- sapply(1:p, function(i) getmin(function(x) wrap_optimDiag(x,i),maxi))
        K2diagonal <-
          sapply(1:p, function(i)
            dichotomie(mini, maxi, function(x)
              wrap_optimDiag(x, i), 1e-5))
        K2diagonalH <-
          sapply(1:r, function(i)
            dichotomie(mini, maxi, function(x)
              wrap_optimDiagH(x, i), 1e-5))

        K2diagonal <- c(K2diagonal, K2diagonalH)
        diag(K2) = K2diagonal

        D2 = diag(K2) %*% t(diag(K2))
        K2[1:p, 1:p] <-
          (1 - sqrt(1 + 4 * Sigma[1:p, 1:p] ^ 2 * D2[1:p, 1:p])) / (2 * Sigma[1:p, 1:p])
        K2[1:p, (p + 1):(p + r)] <-
          (-1 + sqrt(1 + 4 * W ^ 2 * t(D2[1:p, (p + 1):(p + r)]))) / (2 * W)
        K2[(p + 1):(p + r), 1:p] <- t(K2[1:p, (p + 1):(p + r)])
        diag(K2) = diag(K)

        K1 <- nearPD(K2, eig.tol = 1e-6)
        K.tmp = K
        K <- as.matrix(K1$mat)


        Km.tmp <- Km
        if (k == 1) {
          Km <-
            K[1:p, 1:p] - (1 / K[p + 1, p + 1]) * K[1:p, (p + 1):(p + k)] %*% t(K[(p +
                                                                                     1):(p + k), 1:p])
        } else {
          Km <-
            K[1:p, 1:p] - K[1:p, (p + 1):(p + k)] %*% solve(K[(p + 1):(p + k), (p +
                                                                                  1):(p + k)]) %*% K[(p + 1):(p + k), 1:p]
        }
        it = it + 1
        detKm<-max(0.1,min(1e300,det(Km)))
        #print(detKm)
        likelihoods[it] <- -log(detKm) + matrix.trace(Km %*% Sigma[1:p,1:p])-n*log(2*pi)/2
        error[it] = norm(alpha - alpha.tmp,type="F")
        error_[it] = abs(likelihoods[it] - likelihoods[it - 1])
      }
      Sigma.tmp <- Sigma
      Sigma[(p + 1):(p + r), 1:p] <-
        -solve(K[(p + 1):(p + r), (p + 1):(p + r)]) %*% K[(p + 1):(p + r), 1:p] %*% Sigma[1:p, 1:p]
      Sigma[1:p, (p + 1):(p + r)] <-
        -Sigma[1:p, 1:p] %*% K[1:p, (p + 1):(p + r)] %*% solve(K[(p + 1):(p + r), (p + 1):(p + r)])
      Sigma[(p + 1):(p + r), (p + 1):(p + r)] <-
        solve(K[(p + 1):(p + r), (p + 1):(p + r)]) + solve(K[(p + 1):(p + r), (p + 1):(p + r)]) %*%
        K[(p + 1):(p + r), 1:p] %*% Sigma[1:p, 1:p] %*% K[1:p, (p + 1):(p + r)] %*%
        solve(K[(p + 1):(p + r), (p + 1):(p + r)])
      D <- diag(K) %*% t(diag(K))
      d <- ((D - K ^ 2) / D)
      diag(d) = 1
      d = log(d)
      diag(d) = log(diag(K))
      t <- -2 * K[1:p, 1:p] * Sigma[1:p, 1:p]
      diag(t) = -diag(K[1:p, 1:p]) * diag(Sigma[1:p, 1:p])
      f <- computeF(K, Sigma[1:p, 1:p], p)
      m <- matrix(0, p + r, p + r)
      m[1:p, 1:p] <- t
      m[1:p, (p + 1):(p + r)] <- f
      m[(p + 1):(p + r), 1:p] <- t(f)
      gamma <- d + m + log(matrix(rep(pii, (p + r) * (p + r)), p + r))/n
      diag(gamma) = 0
      alpha.tmp = alpha
      magic<-1/n
      alpha <- abs(edge.prob(magic*n*gamma,account.prior = TRUE, q0 = pii))
      #while (sum(is.na(alpha))>0){
      #  magic<-magic/2
      #  alpha <- abs(edge.prob(magic*n*gamma,account.prior = TRUE, q0 = pii))
      #}
      alpha[alpha>1]=1
      alpha[is.na(alpha)]=0
      W <- solve(K[(p + 1):(p + r),(p + 1):(p + r)])%*%K[(p + 1):(p + r), 1:p] %*% S
      V = solve(K[(p + 1):(p + r), (p + 1):(p + r)]) %*% K[(p + 1):(p + r), 1:p]
      V = V %*% Sigma[1:p, 1:p] %*% t(V)
      B <- solve(K[(p + 1):(p + r), (p + 1):(p + r)]) + V
    }
    if (it<max.iter) {
      print("Convergence reached")
    } else{
      print("Maximum number of iterations reached")
    }

    if(r==0){
      P = matrix(0, p, p)
      P = -2 * K* S
      diag(P)= -diag(S * K)
    } else{
      P = matrix(0, p + r, p + r)
      P[1:p, 1:p] = -2 * K[1:p, 1:p] * S
      diag(P)[1:p] = -diag(S * K[1:p, 1:p])
      diag(P)[(p + 1):(p + r)] = -diag(K[(p + 1):(p + r), (p + 1):(p + r)] *
                                         B)
      P[1:p, (p + 1):(p + r)] = (2 * K[1:p, (p + 1):(p + r)] * t(W))
      P[(p + 1):(p + r), 1:p] = t(P[1:p, (p + 1):(p + r)])
    }
    A <- order(abs(alpha) * lower.tri(alpha), decreasing = TRUE)
    A <- A[1:((p + 1) * p / 2)]
    ordered_couples <- lapply(A, function(x)
      returnCouple(x, p + 1))
    return(structure(
      list(
        alpha = alpha,
        ordered_couples = ordered_couples,
        Sigma = Sigma,
        K = K,
        likelihoods = likelihoods,
        gamma = gamma,
        P = P,
        error = error
      )
    ))
  }





tree.EM <- function(S, k, K0, Sigma0, max.iter = 10) {
  # -----------------------------------------------------------------------------------------------------------------------------
  # FUNCTION
  #   Perform EM when k variable are missing and the complete graph is a tree
  # INPUT
  #   S        : emrical covariance matrix of observed variables (pxp matrix)
  #   k        : number of hidden variables (integer)
  #   Sigma0   : initial value of complete covariance matrix ((p+k)x(p+k) matrix)
  #   K0       : initial value complete precision matrix ((p+k)x(p+k) matrix)
  #   max.iter : number of iterations (integer, default is 10)
  # OUTPUT
  #   adjmat    : adjacency matrix describing the tree dependencies ((p+k)x(p+k) matrix)
  #   likelihood: vector containing marginal likelihood at every iteration
  # -----------------------------------------------------------------------------------------------------------------------------
  p <- ncol(S)
  K <- K0
  Tree <- (K != 0) * 1
  Sigma <- Sigma0
  iter <- 1
  if (k == 1) {
    Km <-
      K[1:p, 1:p] - (1 / K[p + 1, p + 1]) * K[1:p, (p + 1):(p + k)] %*% t(K[(p +
                                                                               1):(p + k), 1:p])
  } else {
    Km <-
      K[1:p, 1:p] - K[1:p, (p + 1):(p + k)] %*% solve(K[(p + 1):(p + k), (p +
                                                                            1):(p + k)]) %*% K[(p + 1):(p + k), 1:p]
  }
  likelihood <- rep(0, max.iter)
  while (iter <= max.iter) {
    #likelihood[iter] <- log(det(Km)) - matrix.trace(Km %*% S)
    print(iter)
    ## E step:
    K.tmp <- K
    Tree.tmp <- Tree
    Sigma.tmp <- Sigma
    Sigma[(p + 1):(p + k), 1:p] <-
      -K[(p + 1):(p + k), 1:p] %*% Sigma[1:p, 1:p]
    Sigma[1:p, (p + 1):(p + k)] <-
      -Sigma[1:p, 1:p] %*% K[1:p, (p + 1):(p + k)]
    Sigma[(p + 1):(p + k), (p + 1):(p + k)] <- diag(k)
    + K[(p + 1):(p + k), 1:p] %*% Sigma[1:p, 1:p] %*% K[1:p, (p + 1):(p + k)]

    ## M step:
    Tree <- ChowLiu2(Sigma) #Compute maximum likelihood tree
    E <- which(Tree != 0, arr.ind = T)
    degrees <- rowSums(Tree)
    K.fill <- matrix(0, p + k, p + k)
    for (c in 1:nrow(E)) {
      i <- E[c, 1]
      j <- E[c, 2]
      K.fill[c(i, j), c(i, j)] <-
        K.fill[c(i, j), c(i, j)] + solve(Sigma[c(i, j), c(i, j)])
    }
    K <- K.fill / 2
    diag(K) <- diag(K) - (degrees - 1) / diag(Sigma)
    if (k == 1) {
      Km <- K[1:p, 1:p] - K[1:p, (p + 1):(p + k)] %*% t(K[(p + 1):(p + k), 1:p])
    } else {
      Km <- K[1:p, 1:p] - K[1:p, (p + 1):(p + k)] %*% K[(p + 1):(p + k), 1:p]
    }
    iter <- iter + 1
  }
  adjmat <- (K != 0) * 1
  diag(adjmat) <- 0
  return(structure(list(
    adjmat = adjmat,
    K = K,
    likelihood = likelihood
  )))
}




convex.EM <-
  function(S,
           k,
           lambda,
           K0,
           Sigma0,
           max.iter = 100,
           eps = 1e-3) {
    # -------------------------------------------------------------------------------------------------------------------------------
    # FUNCTION
    #   Perform EM when k variables are missing and perform glasso in M step
    # INPUT
    #   S : empirical covariance matrix of observed variables (pxp matrix)
    #   k : number of hidden variables (integer)
    #   lambda : regularization parameter for graphical lasso in M step (real number)
    #   K0 : initial value of complete precision matrix ((p+k)x(p+k) matrix)
    #   Sigma0 : initial value of complete covariance matrix ((p+k)x(p+k) matrix)
    #   max.iter : number of iterations (integer, default is 10)
    # OUTPUT
    #   adjmat    : adjacency matrix describing the conditional dependencies ((p+k)x(p+k) matrix)
    #   K         : estimator of complete precision matrix ((p+k)x(p+k) matrix)
    #   likelihood: vector containing marginal likelihood at every iteration
    # -------------------------------------------------------------------------------------------------------------------------------
    p <- ncol(S)
    K <- K0
    adjmat <- (K != 0) * 1
    Sigma <- Sigma0
    iter <- 1
    if (k == 1) {
      Km <-
        K[1:p, 1:p] - (1 / K[p + 1, p + 1]) * K[1:p, (p + 1):(p + k)] %*% t(K[(p +
                                                                                 1):(p + k), 1:p])
    } else {
      Km <-
        K[1:p, 1:p] - K[1:p, (p + 1):(p + k)] %*% solve(K[(p + 1):(p + k), (p +
                                                                              1):(p + k)]) %*% K[(p + 1):(p + k), 1:p]
    }
    iter = 1
    error = rep(0, max.iter + 1)
    #    likelihoods <- rep(0, max.iter + 1)
    #    likelihoods[iter] = -log(det(Km)) - matrix.trace(Km %*% S)
    error[iter] = 1
    while (iter <= max.iter && error[iter] > eps) {
      print(paste('EM iteration', iter, sep = " "))

      if (iter == max.iter) {
        print(paste(
          "maximum number of EM iterations reached (",
          max.iter,
          ")",
          sep = ''
        ))
      }
      ## E step:
      K.tmp <- K
      Sigma.tmp <- Sigma
      Sigma[(p + 1):(p + k), 1:p] <-
        -K[(p + 1):(p + k), 1:p] %*% Sigma[1:p, 1:p]
      Sigma[1:p, (p + 1):(p + k)] <-
        -Sigma[1:p, 1:p] %*% K[1:p, (p + 1):(p + k)]
      Sigma[(p + 1):(p + k), (p + 1):(p + k)] <-
        diag(k) + K[(p + 1):(p + k), 1:p] %*% Sigma[1:p, 1:p] %*% K[1:p, (p + 1):(p + k)]

      ## M step:
      mat <- matrix(1, p + k, p + k)
      mat[(p + 1):(p + k), (p + 1):(p + k)] <- 0
      D <- matrix(1, p + k, p + k)
      D[(p + 1):(p + k), (p + 1):(p + k)] <- diag(k)
      lambda_mat <- matrix(lambda, p + k, p + k)

      # Penalization of the whole matrix or not
      #-----------------------------------------
      lambda_mat[(p + 1):(p + k), 1:(p + k)] <- lambda/10
      lambda_mat[1:(p + k), (p + 1):(p + k)] <- lambda/10
      K <-
        glasso(
          Sigma,
          lambda_mat,
          trace = TRUE,
          approx = TRUE,
          penalize.diagonal = FALSE,
          maxit = 500  # Crucial (default Value is 10,000)
        )$wi
      #if (det(K)<0) {stop("Negative Det")}
      if (k > 1) {
        Km <-
          K[1:p, 1:p] - K[1:p, (p + 1):(p + k)] %*% solve(K[(p + 1):(p + k), (p + 1):(p + k)]) %*%
          K[(p + 1):(p + k), 1:p]
      } else{
        Km <- K[1:p, 1:p] - (1 / K[p + 1, p + 1]) * K[1:p, p + 1] %*% t(K[1:p, p +
                                                                            1])
      }
      iter <- iter + 1
      error[iter] = norm(K - K.tmp)
    }
    if (error[iter] <= eps) {
      print("Convergence reached")
    } else
      print("Maximum number of iterations reached")
    mat <- matrix(1, p + k, p + k)
    diag(mat) <- 0
    adjmat <- (K != 0) * 1
    diag(adjmat) <- 0
    return(structure(list(
      adjmat = adjmat,
      K = K    )))
  }


initEM <- function(X = NULL,
#                   S = NULL,
                   cliquelist) {
  # -----------------------------------------------------------------------------------------------------------------------------
  # FUNCTION
  #   Initialize EM by selecting a possible clique and taking the principal component of this clique
  #   as initial value for the hidden variable
  # INPUT
  #   X      : data matrix (nxp)
  #   clique : vector containing the indices of the nodes in the clique. Default is NULL and the clique
  #            is determined with hierarchical clustering.
  # OUTPUT
  #   Sigma0    : initial value of complete covariance matrix ((p+1)x(p+1) matrix)
  #   K0        : initial value of complete precision matrix ((p+1)x(p+1) matrix)
  #   clique    : vector containing the indices of the nodes in the clique
  # -----------------------------------------------------------------------------------------------------------------------------
  n<-nrow(X)
  num_clusts <- length(cliquelist)
  pr <-
    lapply(1:num_clusts, function(k)
      prcomp(t(X[, cliquelist[[k]]]), center = TRUE, scale = TRUE))
  prmat <-
    lapply(1:num_clusts, function(k)
      scale(pr[[k]]$rotation[, 1], center = TRUE, scale = TRUE))
  Y <- scale(X, center = TRUE, scale = TRUE)
  Y <-
    scale(cbind(Y, do.call(cbind, prmat)), center = TRUE, scale = TRUE)
  mu <- matrix(0, ncol(Y), 1)
  Y <- Y + mvrnorm(n, mu, 0.01 * diag(ncol(Y)))
  Sigma0 <- 1 / n * t(Y) %*% (Y)
  K0 <- pinv(Sigma0)
  return(structure(list(
    Sigma0 = Sigma0,
    K0 = K0,
    cliquelist = cliquelist
  )))
}


findCliques<-function(X,k){
  # -----------------------------------------------------------------------------------------------------------------------------
  # FUNCTION
  # Find k groups of nodes which are 'clique likes"
  # INPUT
  #   X      : data matrix (nxp)
  #   k      : number of groups
  #            is determined with hierarchical clustering.
  # OUTPUT
  #   cliquelist    : a list of groups of indices
  # -----------------------------------------------------------------------------------------------------------------------------

  Y = 1- (cor(X))^2
  res.kmeans<-kmeans(Y,k)
  # Compute the mean distance within each cluster (clique should have the smallest mean distances)
  scores<-res.kmeans$withinss / ((res.kmeans$size)*(res.kmeans$size-1)/2)
  # Renumbering of the cluster according the score (cluster one will have the smallest score)
  reordered.cluster<-order(scores)[res.kmeans$cluster]
  return(cliquelist = split(1:length(reordered.cluster),reordered.cluster))
}

matrix.trace<-function(M){sum(diag(M))}


get.lambda.l1 <- function(S) {
  r.max <- max(0,max(abs(S-diag(diag(S)))))
  return(r.max)
}



removeDiag<-function(M){
  diag(M)<-NA
  return(matrix(M[!is.na(M)],nrow(M),ncol(M)-1))
}

ChowLiu2 <- function(Sigma){
  # ---------------------------------------------------------------------------------------------------------
  # FUNCTION
  #   Build Chow-Liu tree from empirical covariance matrix
  # INPUT
  #   Sigma    : empirical covariance matrix (pxp matrix)
  # OUTPUT
  #   adjmat   : adjacency matrix describing conditional dependence structure
  # ---------------------------------------------------------------------------------------------------------
  df <- data.frame(Sigma)
  colnames(df) <- c(1:nrow(Sigma))
  graph <- chow.liu(df)
  adjmat <- amat(graph)
  return(adjmat)
}


pinv<-function (A, tol = .Machine$double.eps^(2/3))
{
  stopifnot(is.numeric(A), length(dim(A)) == 2, is.matrix(A))
  s <- svd(A)
  p <- (s$d > max(tol * s$d[1], 0))
  if (all(p)) {
    mp <- s$v %*% (1/s$d * t(s$u))
  }
  else if (any(p)) {
    mp <- s$v[, p, drop = FALSE] %*% (1/s$d[p] * t(s$u[,
                                                       p, drop = FALSE]))
  }
  else {
    mp <- matrix(0, nrow = ncol(A), ncol = nrow(A))
  }
  return(mp)
}



returnCouple <- function(ind, p){
  # ---------------------------------------------------------------------------------------------------------
  # FUNCTION
  #   Convert linear matrix index into row/column indices
  # INPUT
  #   ind    : matrix index (integer in 1:p^2
  #   p      : size of matrix (integer)
  # OUTPUT
  #   couple : vector of size 2 containing the row and column indices to get the ind_th value of a matrix of
  #            size p
  # ---------------------------------------------------------------------------------------------------------

  if(ind%%p==0){
    i <- ind%/%p
    j <- p
  } else{
    i <- ind%/%p+1
    j <- ind%%p
  }
  couple <- c(i,j)
  return(couple)
}



computeF <- function(K, S, p){
  # ---------------------------------------------------------------------------------------------------------
  # FUNCTION
  #   Compute the log(f_ih), i€O, h€H coefficients of the posterior edges probabilities
  # INPUT
  #   K        : (p+rxp+r matrix) p number of observed variables r number of hidden variables
  #   n        : number of samples
  #   S        : Empirical covariance matrix of observed variables (pxp matrix)
  # OUTPUT
  #   f        : pxr matrix containing the coefficients log(f_ih). f[i,h]=log(f_ih)
  # ---------------------------------------------------------------------------------------------------------

  if(is.null(ncol(K))){
    r <- length(K)-p
  } else{
    r <- ncol(K)-p
  }
  f <- (K[1:p,(p+1):(p+r)]%*%solve(K[(p+1):(p+r),(p+1):(p+r)]))*t(t(K[1:p,(p+1):(p+r)])%*%S)
  f <- f
  return(f)
}

optimDiag <- function(i,k_ii, K, Sigma, alpha){
  p <- nrow(K)
  k <- K[i,setdiff(1:p,i)]
  a <- alpha[i,setdiff(1:p,i)]
  return(1/k_ii-Sigma[i,i]+sum(a*k^2/(k_ii*diag(K)[setdiff(1:p,i)]-k^2)))
}


dichotomie <- function(a, b, f, epsilon){
  # ---------------------------------------------------------------------------------------------------------
  # FUNCTION
  #   find the zero of monotonous function f by binary search.
  # INPUT
  #   a        : minimum argument of f
  #   b        : maximum argument of f
  #   epsilon  : tolerance (length of last interval)
  # OUTPUT
  #   c        : approximate zero of f with epsilon tolerance
  # ---------------------------------------------------------------------------------------------------------

  min <- a
  max <- b
  sgn <- sign(f(max)-f(min))
  c <- (max+min)/2
  while(abs(f(c))>1e-5){
    c <- (max+min)/2
    if(sgn*f(c)>0){
      max <- c
    } else{
      min <- c
    }
  }
  return(c)
}

getmin <- function(f,max){
  vals=f(seq(0,1,length.out=1000))
  mini=min(intersect(which(abs(vals)<Inf), which(sign(vals)+sign(f(max))==0)))*1e-3
}



