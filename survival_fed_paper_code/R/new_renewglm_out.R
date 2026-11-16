#X = X_mat; y = y; type = "binomial"; betahat = betahat; infomats = infomats; intercept = FALSE; s = NULL; phi = NULL
RenewGLM_out_ori <- function(X, y, type, betahat, infomats, intercept, s, phi){
    tol <- 1e-6;
    max_iter <- 1000;
    
    p <- dim(as.matrix(X))[2];
    
    betahat_old <- betahat;
    
    W <- invlinkdiv(X, betahat_old, type=type)
    H <- cp(X,y,W)
    
    #eps <- 1e-6
    #matrix <- infomats + H + diag(eps, p)
    matrix <- infomats + H 
    
    U=chol(matrix)
    L=t(U)
    
    for (r in 1:max_iter){
      
      g_0=t(crossprod((y-invlink(X, betahat, type)),X));
      g_1=-t(crossprod((betahat-betahat_old),infomats));
      g <- g_0 + g_1;
      
      d_beta <- backsolve(U,forwardsolve(L,g))
      
      df_beta <- crossprod(g,d_beta);
      if (abs(df_beta)<tol){
        break
      }else {
        betahat <- betahat + drop(d_beta);
      }
      }
    # if(type=="gaussian"){
    #   s <- s+length(y)
    #   phi<-(s-length(y)-p)/(s-p)*phi+1/(s-p)*t(betahat_old)%*%infomats%*%(betahat_old-betahat)+t(y)%*%(y-X%*%betahat)/(s-p);
    # }
    
    W <- invlinkdiv(X,betahat,type=type)
    H_new <- cp(X,y,W)
    infomats<-infomats+H_new
    
    # if(type=="gaussian"){
    #   return(list(betahat,infomats,s,phi))
    # }else{
    #   return(list(betahat,infomats))
    # }
    return(list(betahat,infomats))
  }

RenewGLM_out <- function(X, y, id, type, betahat, infomats, meat_list, intercept = FALSE, s = NULL, phi = NULL) {
  tol <- 1e-6
  max_iter <- 1000
  
  p <- dim(as.matrix(X))[2]
  
  betahat_old <- betahat
  
  # Parameter update
  W <- invlinkdiv(X, betahat_old, type = type)
  H <- cp(X, y, W)
  
  eps <- 1e-6
  matrix <- infomats + H + diag(eps, p)
  
  U <- chol(matrix)
  L <- t(U)
  
  for (r in 1:max_iter) {
    print(r)
    g_0 <- t(crossprod((y - invlink(X, betahat, type)), X))
    g_1 <- -t(crossprod((betahat - betahat_old), infomats))
    g <- g_0 + g_1
    
    d_beta <- backsolve(U, forwardsolve(L, g))
    df_beta <- crossprod(g, d_beta)
    
    if (abs(df_beta) < tol) {
      break
    } else {
      betahat <- betahat + drop(d_beta)
    }
  }
  
  W <- invlinkdiv(X, betahat, type = type)
  H_new <- cp(X, y, W)
  infomats <- infomats + H_new
  
  # Robust variance computation (sandwich estimator)
  unique_ids <- unique(id)
  Meat <- matrix(0, p, p)
  
  for (cluster in unique_ids) {
    idx <- which(id == cluster)
    X_cluster <- X[idx, , drop = FALSE]
    y_cluster <- y[idx]
    
    mu <- invlink(X_cluster, betahat, type)
    residuals <- y_cluster - mu
    score_cluster <- t(X_cluster) %*% residuals
   
    Meat <-  Meat + (score_cluster %*% t(score_cluster)) 
  }
  
  return(list(betahat = betahat, infomats = infomats, meat_cum = Meat))
}

cp <-
  function(X, y, w){
    if (length(y)==1){
      H <- w*tcrossprod(X,X)
    }else{
      H <- crossprod(sqrt(w)*X)
    }
  }

invlink <-
  function(X, beta, type){
    #  if (length(beta)==1){
    #    eta <- X * beta
    #  }else{
    eta <- drop(X %*% beta)
    #  }
    if(type=="gaussian"){ out <- eta }
    if(type=="binomial"){ out <- exp(eta)/(1 + exp(eta)) }
    if(type=="poisson") { out <- exp(eta) }
    if(type=="cloglog"){out <- 1 - exp(-exp(eta))}
    out
  }

invlinkdiv <-
  function(X, beta, type){
    eta <- drop(X %*% beta)
    if(type=="gaussian"){ out <- 1 }
    if(type=="binomial"){ out <- exp(eta)/(1 + exp(eta))^2 }
    if(type=="poisson") { out <- exp(eta) }
    if(type=="cloglog"){out <- exp(eta-exp(eta))}
    out
}



