#library(saturnin)
critere <- function(X, EM.res,pii=0.5){
  K=EM.res$K
  Sigma=EM.res$Sigma
  log_gamma=EM.res$gamma
  alpha=EM.res$alpha
  P=EM.res$P
  diag(log_gamma)=0
  n=nrow(X)
  q=nrow(log_gamma)
  p=ncol(X)
  r=nrow(Sigma)-p
  #pii=log(matrix(2/(p+r),p+r,p+r))/n
  pii_mat=log(matrix(pii,p+r,p+r))
  diag(pii_mat)=0
  diag(P)=0
  delta=-exp(log_gamma)
  diag(delta)=rowSums(exp(log_gamma))
  H.T=log(det(delta[2:q,2:q]))-sum(sum(log_gamma*alpha)) # Entropy of trees
  #H.T=sum(log(d))-sum(sum(log_gamma*alpha)) # Entropy of trees
  if (r==1){
  	E.H.Xh=r*log(2*pi*exp(1))/2-(1/2)*sum(log(K[(p+1),(p+1)])) # Entropy of hidden variables
  } else E.H.Xh=r*log(2*pi*exp(1))/2-(1/2)*sum(log(diag(K[(p+1):(p+r),(p+1):(p+r)]))) # Entropy of hidden variables
  pen.r=(log(n))*(p*(p+1)/2+r*p+r)
  delta0 = -exp(pii_mat)+diag(rowSums(exp(pii_mat)))
  Z0 = det(delta0[2:q,2:q])
  #mat=alpha*(pii-Z0+P- (p+r)*log(2*pi)/2)
  mat=alpha*(pii_mat+P)
  E.c=sum(sum(mat)) - Z0-(p+r)*log(2*pi)/2# Expectation of complete likelihood
  logPY = E.c+E.H.Xh+H.T # Marginal log-likelihood of observed variables
  return(list(ICL_T=logPY-H.T-pen.r, ICL_ZT=E.c-pen.r, BIC=logPY-pen.r))
}
