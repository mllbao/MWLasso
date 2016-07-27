choose_X <- function(G, n) {
  if (n != nrow(G)){
    row = sort(sample(x=1:nrow(G), size = n, replace = F))
    X = G[row,]; 
  }else if(n == nrow(G)){
    X = G;
  }
  sdg = rep(0,ncol(X));
  for ( i in 1:ncol(X) ){
    g = X[,i]
    sdg[i] = sd(g, na.rm=T)
  }
  p = ncol(X);
  index_1 = (1:p)[sdg!=0]; index_0 = (1:p)[sdg==0]
  X_1 = X[, index_1]; X_0 = X[, index_0];
  result = list("X_1"=X_1, "X_0"=X_0, "index_1"=index_1, "index_0"=index_0);
  return (result);
}

get_Nk <- function(List){
  X_1 = List$X_1; index_1 = List$index_1;
  X_0 = List$X_0; index_0 = List$index_0;
  p = ncol(X_1); q = ncol(X_0);
  X = matrix(0, nrow = nrow(X_1), ncol = p + q); 
  X[, index_1] = X_1;
  if (q!=0){
    X[, index_0] = X_0;
  }
  n = ncol(X);
  Nk = rep(0, n);
  for ( i in 1:n ){
    g = X[,i]
    Nk[i] = length(g[!is.na(g)]);
  }
  return (Nk);
}

Stand_X <- function (List) {
  X_1 = List$X_1; index_1 = List$index_1;
  X_0 = List$X_0; index_0 = List$index_0;
  for ( i in 1:ncol(X_1) ){
    g = X_1[,i]; n = length(X_1[,i][!is.na(X_1[,i])]);
    X_1[,i] = sqrt(n/(n-1))*(g - mean(g, na.rm=T))/sd(g, na.rm=T)
  }  
  if (ncol(X_0)!=0){
    for ( i in 1:ncol(X_0) ){
      X_0[,i] = rep(0, length(X_0[,i]))
    } 
  }
  List$X_1 = X_1; List$X_0 = X_0;
  return (List);
}

cor_coef <- function (d, List) {
  X = List$X_1; p = ncol(X);
  W = matrix(0, nrow = d, ncol = p-d+1);
  for ( i in 1:d ){
    W[i,] = c(i:(i+p-d))
  }
  coef = NULL; V = NULL;
  for ( k in 1:p ) {  
    W_k = W[, seq(max(1, k-d+1), min(k, p-d+1), 1)];
    
    V_k = as.vector(W_k);
    V_k = V_k[V_k!=k];
    quasi = rep(0, length(V_k));
    for ( i in 1:length(quasi) ) {
      q = sum(X[,V_k[i]]*X[,k], na.rm = T)/sqrt(sum(X[,V_k[i]]*X[,V_k[i]], na.rm = T)*sum(X[,k]*X[,k], na.rm = T))
      quasi[i] = abs(q)
    }
    coef = c(coef, list(quasi))
    V = c(V, list(V_k))
  }
  result = list(coef = coef, V = V)
  return (result)
}

CDM <- function (List, Y, Nk, cc, lambda, eta, d, epson = 10^(-10), M = 100) {
  X_1 = List$X_1; index_1 = List$index_1;
  X_0 = List$X_0; index_0 = List$index_0;
  p = ncol(X_1); q = ncol(X_0);
  beta = rep(0,p);
  P = rep(0, p); Q = rep(0, p); S = rep(0, p);
  coef = cc$coef; V = cc$V; Nn = Nk[index_1];
  for ( k in 1:p ) {  
    quasi = coef[[k]];
    P[k] = 1/2*(1/Nn[k]*sum(na.omit(X_1[,k]^2)) + eta*sum(quasi))
    Q[k] = -1/Nn[k] * sum(na.omit(Y*X_1[,k]))
  }
  err = 1; ite = 1;
  while ( (err >= epson) && (ite <= M) ) {
    temp = beta; 
    for ( k in 1:p ) {  
      quasi = coef[[k]];
      beta_k = abs(beta[V[[k]]])
      S[k] = lambda - eta*sum(quasi*beta_k)
      beta[k] = -sign(Q[k])*(max(abs(Q[k])-S[k], 0))/(2*P[k])
    }
    ite = ite + 1; 
    err = sum(abs(temp - beta)); 
  }
  result = rep(0, p+q);
  result[index_1] = beta;
  return (result);
}

pred_num <- function(List, Y, Nk, cc, d, lambda, eta, epson = 10^(-10), M = 100) {
  eta = eta/(d-1);
  beta = CDM(List=List, Y=Y, Nk=Nk, cc=cc, lambda=lambda, eta=eta, d=d, epson = epson, M = M);
  result = length(beta[beta!=0]);
  return (result)
}

find_lambda <- function(List, Y, Nk, cc, d, number, gamma2, epson = 10^(-10), M = 100){
  X_1 = List$X_1; index_1 = List$index_1;
  p = ncol(X_1)
  max_lambda = rep(0, p)
  for ( i in 1:p ) {
    max_lambda[i] = sum(na.omit(X_1[,i]*Y))/Nk[i]
  }
  lambda_ub = max(max_lambda) 
  max_gamma1 = lambda_ub/gamma2
  min_gamma1 = 0.1*max_gamma1
  gamma1_l = min_gamma1; gamma1_u = max_gamma1;
  lambda_l = gamma1_l*gamma2; eta_l = gamma1_l-lambda_l;
  lambda_u = gamma1_u*gamma2; eta_u = gamma1_u-lambda_u;
  tau_l = pred_num(List, Y, Nk, cc, d, lambda_l, eta_l, epson, M);
  tau_u = pred_num(List, Y, Nk, cc, d, lambda_u, eta_u, epson, M);
  if(tau_l == number){
    result = lambda_l
  }else if(tau_u == number){
    result = lambda_u
  }
  k = 1;
  while((tau_l > number)&&(tau_u < number)&&(k <= 50)){
    gamma1_m = 0.5*(gamma1_l+gamma1_u);
    lambda_m = gamma1_m*gamma2; eta_m = gamma1_m-lambda_m;
    tau_m = pred_num(List, Y, Nk, cc, d, lambda_m, eta_m, epson, M);
    if(tau_m < number){
      gamma1_u = gamma1_m
    }else if(tau_m > number){
      gamma1_l = gamma1_m
    }else if(tau_m == number){
      result = lambda_m; break;
    }
    k = k+1;
    if (k == 50)
      result = lambda_m;
  }
  return (result)
}

get_beta0 <- function(List, Y) {
  X_1 = List$X_1; index_1 = List$index_1;
  X_0 = List$X_0; index_0 = List$index_0;
  X = X_1;
  p = ncol(X);
  beta0 = rep(0, p);
  for ( k in 1:p ){
    x = X[,k];
    y = Y[!is.na(x)]; x = x[!is.na(x)];
    model = glm(formula = y ~ x, family = binomial(logit));
    beta0[k] = model$coefficients[1];
  }
  return (beta0);
}

CDM_qualitative <- function(List, Y, Nk, cc, initial.beta0, lambda, eta, d, epson = 10^(-8), M = 100){
  X_1 = List$X_1; index_1 = List$index_1;
  X_0 = List$X_0; index_0 = List$index_0;
  X = X_1;
  p = ncol(X); n = nrow(X); q = ncol(X_0);
  P = rep(0, p); Q = rep(0, p); S = rep(0, p); 
  Prob = matrix(0, ncol = p, nrow = n); Wei = Prob; Z = Prob;
  coef = cc$coef; V = cc$V; Nn = Nk[index_1];
  beta0 = initial.beta0;
  beta1 = rep(0, p);
  err = 1; ite = 1;
  while ((err >= epson) && (ite <= M)) {
    temp = Prob; 
    for ( j in 1:p ) {
      x = X[,j];
      y = Y[!is.na(x)]; x = x[!is.na(x)];
      prob = exp(beta0[j] + beta1[j]*x)/(1 + exp(beta0[j] + beta1[j]*x)); Prob[,j][!is.na(X[,j])] = prob;
      z = beta0[j] + beta1[j]*x + (y - prob)/(prob * (1 - prob)); Z[,j][!is.na(X[,j])] = z;
      wei = prob*(1-prob); Wei[,j][!is.na(X[,j])] = wei;
    }
    for ( j in 1:p ){
      x = X[,j];
      quasi = coef[[j]];
      P[j] = 1/2*(1/Nn[j]*sum(na.omit(Wei[,j] * x^2)) + eta*sum(quasi))
      Q[j] = -1/Nn[j] * sum(na.omit(Wei[,j] * (Z[,j] - beta0[j]) * x))
    }  
    err.in = 1; ite.in = 1;
    while ( (err.in >= epson) && (ite.in <= M) ) {
      temp.in = beta1; 
      for ( k in 1:p ) { 
        quasi = coef[[k]];
        beta1_k = abs(beta1[V[[k]]])
        S[k] = lambda - eta*sum(quasi*beta1_k)
        beta1[k] = -sign(Q[k])*(max(abs(Q[k])-S[k], 0))/(2*P[k])
      }
      ite.in = ite.in + 1; 
      err.in = max(abs(temp.in - beta1)); 
    }
    ite = ite + 1; 
    err = max(abs(temp - Prob));
  }
  result = rep(0, p+q);
  result[index_1] = beta1;
  return (result);
}


pred_num_qualitative <- function(List, Y, Nk, cc, initial.beta0, d, lambda, eta, epson = 10^(-8), M = 100) {
  eta = eta/(d-1);
  beta = CDM_qualitative(List=List, Y=Y, Nk=Nk, cc=cc, initial.beta0=initial.beta0, lambda=lambda, eta=eta, d=d, epson = epson, M = M);
  result = length(beta[beta!=0]);
  return (result)
}

find_lambda_qualitative <- function(List, Y, Nk, cc, initial.beta0, d, number, gamma2, epson = 10^(-8), M = 100){
  X_1 = List$X_1; index_1 = List$index_1;
  p = ncol(X_1);
  max_lambda = rep(0, p)
  for ( i in 1:p ) {
    max_lambda[i] = sum(na.omit(X_1[,i]*Y))/Nk[i]
  }
  lambda_ub = max(max_lambda) 
  max_gamma1 = lambda_ub/gamma2
  min_gamma1 = 0.1*max_gamma1
  gamma1_l = min_gamma1; gamma1_u = max_gamma1;
  lambda_l = gamma1_l*gamma2; eta_l = gamma1_l-lambda_l;
  lambda_u = gamma1_u*gamma2; eta_u = gamma1_u-lambda_u;
  tau_l = pred_num_qualitative(List, Y, Nk, cc, initial.beta0, d, lambda_l, eta_l, epson, M);
  tau_u = pred_num_qualitative(List, Y, Nk, cc, initial.beta0, d, lambda_u, eta_u, epson, M);
  if(tau_l == number){
    result = lambda_l
  }else if(tau_u == number){
    result = lambda_u
  }
  while((tau_l > number)&&(tau_u < number)){
    gamma1_m = 0.5*(gamma1_l+gamma1_u);
    lambda_m = gamma1_m*gamma2; eta_m = gamma1_m-lambda_m;
    tau_m = pred_num_qualitative(List, Y, Nk, cc, initial.beta0, d, lambda_m, eta_m, epson, M);
    if(tau_m < number){
      gamma1_u = gamma1_m
    }else if(tau_m > number){
      gamma1_l = gamma1_m
    }else if(tau_m == number){
      result = lambda_m; break;
    }
  }
  return (result)
}

MWLasso <- function(X, Y, lambda, eta, d, method='linear', epson = 10^(-10), M = 100){
  n = nrow(X)
  List = choose_X(X, n)
  List = Stand_X(List)
  Nk = get_Nk(List)
  cc = cor_coef(d, List)
  eta = eta/(d-1)
  if (method == 'linear'){
    Y = Y - mean(Y)
    beta.hat = CDM(List, Y, Nk, cc, lambda, eta, d, epson, M)
    return (beta.hat)
  } else if(method == 'logistic'){
    initial.beta0 = get_beta0(List, Y)
    beta.hat = CDM_qualitative(List, Y, Nk, cc, initial.beta0, lambda, eta, d, epson, M)
    return (beta.hat)
  } else {
    print("Error: method should be linear or logistic!")
  }
}

MW_parameters <- function(X, Y, d, number, gamma2, method="linear", epson = 10^(-10), M = 100){
  n = nrow(X)
  List = choose_X(X, n)
  List = Stand_X(List)
  Nk = get_Nk(List)
  cc = cor_coef(d, List)
  if (method == 'linear'){
    Y = Y - mean(Y)
    lambda = find_lambda(List, Y, Nk, cc, d, number, gamma2, epson, M)
    eta = lambda*(1/gamma2 - 1)
    cat('Tuning parameters are: lambda = ',lambda,', eta = ', eta, '.', sep = "")
    return (c(lambda, eta))
  } else if (method == 'logistic'){
    initial.beta0 = get_beta0(List, Y)
    lambda = find_lambda_qualitative(List, Y, Nk, cc, initial.beta0, d, number, gamma2, epson, M)
    eta = lambda*(1/gamma2 - 1)
    cat('Tuning parameters are: lambda = ',lambda,', eta = ', eta, '.', sep = "")
    return (c(lambda, eta))
  } else {
    print ("Error: method should be linear or logistic!")
  }
}
