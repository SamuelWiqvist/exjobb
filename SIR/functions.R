# IF1_SIR
IF1_SIR <- function(theta_0,Y_hat, M, N,J_seq,tau_seq,a_seq,sigma_seq,Sigma){
  
  # Description
  # The IF2 algorithm based upon the algorithm described in
  # Ionides (2011)
  #
  # Input:
  # theta0 = start value 
  # Y_hat = observatrions 
  # M = nbr of iterations 
  # N = nbr of timestep 
  # J_seq = nbr of particles
  # tau_seq = cooling factor for theta
  # sigma_seq = cooling factor for theta
  # a_seq = cooling factor for the updating formula
  # Sigma = covariance matrix for the distribtuion k 
  #
  # Output:
  # theta_obs = the estimated values of theta for each M-iteration
  # log_lik_vec = estiamtion of the log-lik
    
  # Print msg 
  print('Starting IF1_SIR')
    
  log_lik_vec = rep(0,M)
    
  theta_vec = matrix(0,M+1,3)
  theta_vec[1,] = theta_0
    
  # set distribution for parameter perturbation 
  limit = 1000;
  pd_beta_mu = 0;
  pd_beta_sigma = Sigma[1,1];
  # kappa_beta = rtruncnorm(n,-limit,limit, pd_beta_mu, pd_beta_sigma) 

  pd_gamma_mu = 0;
  pd_gamma_sigma = Sigma[2,2];
  # kappa_beta = rtruncnorm(n,-limit,limit, pd_gamma_mu, pd_gamma_sigma) 
  
  pd_rho_transformed_mu = 0;
  pd_rho_transformed_sigma = Sigma[3,3];
  # kappa_beta = rtruncnorm(n,-limit,limit, pd_rho_transformed_mu, pd_rho_transformed_sigma) 
  
  # pre-allocate vectors 
  theta_bar_vec = matrix(0, N, 3);
  theta_var_vec = matrix(0, N, 3);
  
 
  for(m in 1:M){
    print(m)
    # Set J
    J = J_seq[m];
    
    theta_F = matrix(0,J,3)
    theta_P = matrix(0,J,3)
    
    # set marix for weigths
    
    w_m = matrix(0, N, J)
    
    # set indices
    index  = 1:J;
    
    # draw parameter perturbations
    Z_beta = rtruncnorm(J,-limit,limit, pd_beta_mu, pd_beta_sigma)
    Z_gamma = rtruncnorm(J,-limit,limit, pd_gamma_mu, pd_gamma_sigma) 
    Z_rho = rtruncnorm(J,-limit,limit, pd_rho_transformed_mu, pd_rho_transformed_sigma) 
    
    # Initilize parameters
    tau_Z_beta = tau_seq[m]*Z_beta
    tau_Z_gamma = tau_seq[m]*Z_gamma
    tau_Z_rho = tau_seq[m]*Z_rho
    
    theta_F[,1] = theta_vec[m,1] + tau_Z_beta
    theta_F[,2] = theta_vec[m,2] + tau_Z_gamma
    theta_F[,3] = theta_vec[m,3] + tau_Z_rho
    
    # Initilize states
    X_F = round( rmvnorm(J,c(30000, 800, 470000, 400), diag(c(50,30,500,10)) ) ) # X_F = [S I R H]
    
    popsize = rowSums(X_F[,1:3]); 
    
    # calc theta_bar_0
    theta_bar_0 = colMeans(theta_F);
    
    for(n in 1:N){
      # Perturb parameters
      Z_beta = rtruncnorm(J,-limit,limit, pd_beta_mu, pd_beta_sigma)
      Z_gamma = rtruncnorm(J,-limit,limit, pd_gamma_mu, pd_gamma_sigma) 
      Z_rho = rtruncnorm(J,-limit,limit, pd_rho_transformed_mu, pd_rho_transformed_sigma) 
      
      theta_P[,1] = theta_F[,1] + sigma_seq[m]*Z_beta
      theta_P[,2] = theta_F[,2] + sigma_seq[m]*Z_gamma
      theta_P[,3] = theta_F[,3] + sigma_seq[m]*Z_rho
      

      # Simulate prediction particles
      X_P = simulatior_X(theta_P,X_F, popsize) # need to fix this

      # Evaluate weigths
      w = evaluator_Y_X(Y_hat[n], X_P, theta_P)
      prob_resample = w/sum(w) # normalize 
      # store wegiths
      
      w_m[n,] = w 
      
      # Systematic resampling of indecies 
      K = sample(index, J, replace = TRUE, prob = prob_resample)

      # Resaple 
      X_F = X_P[K,]
      theta_F = theta_P[K,]
      
      # calc mean
      theta_bar_vec[n, ] = colMeans(theta_F);
      
      if(n == 1){
        theta_var_vec[n,] = (1/(J-1))*c(sum((theta_P[,1] - theta_bar_0[1])^2), 
                                      sum((theta_P[,2] - theta_bar_0[2])^2),
                                      sum((theta_P[,3] - theta_bar_0[3])^2));
      }else{
        theta_var_vec[n,] = (1/(J-1))*c(sum((theta_P[,1] - theta_bar_vec[n-1,1])^2),  
                                        sum((theta_P[,2] - theta_bar_vec[n-1,2])^2), 
                                        sum((theta_P[,3] - theta_bar_vec[n-1,3])^2));
      }
      
    }
    # calc log lik
    log_lik_vec[m] = sum(log(rowMeans(w_m)));
    
    # update 
    theta_var_vec = theta_var_vec^(-1);
    end_val = nrow(theta_bar_vec)
    G_m_theta = c(sum( theta_var_vec[,1]*(theta_bar_vec[1:N,1] - c(theta_bar_0[1], theta_bar_vec[1:(end_val-1),1]))),
                 sum( theta_var_vec[,2]*(theta_bar_vec[1:N,2] - c(theta_bar_0[2], theta_bar_vec[1:(end_val-1),2]))),
                 sum( theta_var_vec[,3]*(theta_bar_vec[1:N,3] - c(theta_bar_0[3], theta_bar_vec[1:(end_val-1),3]))))
    
    
    theta_vec[m+1,] = theta_vec[m,] + a_seq[m,]*G_m_theta;
    
  
  }
  

  # create output list
  output_list = list('theta_IF1' = theta_vec, 'loglik' = log_lik_vec)
  
  return(output_list)
}

# IF2_SIR
IF2_SIR <- function(theta0,Y_hat, M, sigma_m,J,N){
  # Description
  # Estimates the parameters theta using IF2 agirithm described in 
  # Ionides (2015)
  #
  # Inputs: 
  # theta0 = inital particle swarm
  # y_hat = obseravtions 
  # M = nbr of iterations
  # sigma_m: preturbation sequance (decreses exponentially)
  # J: nbr of particles
  # N: nbr of time steps 
  
  Sigma = diag(3);
  Sigma[1,1] = 10
  Sigma[2,2] = 1
  Sigma[3,3] = 0.01
  
  # print msg 
  
  print('Starting IF2_SIR')
  
  
  # simulator for f_x0 : N(0,1), 
  # simulator for f_xn | x_n-1 ; omega : N(omgea*x_n-1,1) and
  # simulator for f yn | xn : normal(xn,2) are hard coded 
  # h_n: gaussian with mean omega_j and variance sigma_m
  
  # pre-allocate matrix for particels
  theta_m = matrix(0,M,J*3);  
  
  # set indecies
  index = 1:J;  
  
  # pre-allocationg matrices
  theta_F = matrix(0,J,3)
  theta_P = matrix(0,J,3)
  theta_F_minus = matrix(0,J,3)
  # procidure
  for (m in 1:M) { # loop for nbr of iterations 
    print(m)
    # initilize parameters
    if( m == 1 ){ 
      for( i in 1:nrow(theta0)){
        theta_F[i,] = rmvnorm(1,theta0[i,], sigma_m[m]*Sigma)
      }
    }else{
      theta0[,1] = theta_m[m-1,(1:J)]
      theta0[,2] = theta_m[m-1,(J+1):(2*J)]
      theta0[,3] = theta_m[m-1,(2*J+1):(2*J)]
      for( i in 1:nrow(theta0)){
        theta_F[i,] = rmvnorm(1,theta0[i,], sigma_m[m]*Sigma)
      }
    }
                       
    # Initilize states
    X_F = round( rmvnorm(J,c(30000, 800, 470000, 400), diag(c(50,30,500,10)) ) ) # X_F = [S I R H]
    
    popsize = rowSums(X_F[,1:3]); 
    
    for( n in 1:N) { # loop for nbr of times steps 
                     
       # Preturb parameters
       for( i in 1:nrow(theta_F) ){
        theta_P[i,] = rmvnorm(1,theta_F[i,], sigma_m[m]*Sigma); 
       }
       # Simulate prediction particles
       X_P = simulatior_X(theta_P,X_F, popsize) # need to fix this
       
       # Evaluate weigths
       w = evaluator_Y_X(Y_hat[n], X_P, theta_P)
       prob_resample = w/sum(w) # normalize 
       
       # Systematic resampling of indecies
       K = sample(index, J, replace = TRUE, prob = prob_resample)
       
       # Resaple 
       theta_F = theta_P[K,];
       X_F = X_P[K,];
    }
    # Store new parameter swarm
    theta_m[m,] = c(theta_F[,1], theta_F[,2], theta_F[,3]) ;
  }
  
  # Add the inital parameter swarm
  return(theta_m)
}             
             
# AIF_SIR
AIF_SIR <- function(theta_0,Y_hat, M, N,J_seq,tau_seq,sigma_seq,Sigma, a_m, c_m, burn_in){
  
  # Description
  # The IF2 algorithm based upon the algorithm described in
  # Ionides (2011)
  #
  # Input:
  # theta0 = start value 
  # Y_hat = observatrions 
  # M = nbr of iterations 
  # N = nbr of timestep 
  # J_seq = nbr of particles
  # tau_seq = cooling factor for theta
  # sigma_seq = cooling factor for theta
  # a_seq = cooling factor for the updating formula
  # Sigma = covariance matrix for the distribtuion k 
  #
  # Output:
  # theta_obs = the estimated values of theta for each M-iteration
  # log_lik_vec = estiamtion of the log-lik
  
  # Print msg 
  print('Starting AIF_SIR')
  
  log_lik_vec = rep(0,M)
  
  theta_vec = matrix(0,M+1,3)
  theta_vec[1,] = theta_0
  
  # set distribution for parameter perturbation 
  limit = 1000;
  pd_beta_mu = 0;
  pd_beta_sigma = Sigma[1,1];
  # kappa_beta = rtruncnorm(n,-limit,limit, pd_beta_mu, pd_beta_sigma) 
  
  pd_gamma_mu = 0;
  pd_gamma_sigma = Sigma[2,2];
  # kappa_beta = rtruncnorm(n,-limit,limit, pd_gamma_mu, pd_gamma_sigma) 
  
  pd_rho_transformed_mu = 0;
  pd_rho_transformed_sigma = Sigma[3,3];
  # kappa_beta = rtruncnorm(n,-limit,limit, pd_rho_transformed_mu, pd_rho_transformed_sigma) 
  
  # pre-allocate vectors 
  theta_bar_vec = matrix(0, N, 3);
  theta_bar_vec_add = matrix(0, N, 3);
  theta_bar_vec_sub = matrix(0, N, 3);
  
  # vectors to store the variance of theta values in 
  theta_var_vec = matrix(0, N, 3);
  theta_var_vec_add = matrix(0, N, 3);
  theta_var_vec_sub = matrix(0, N, 3);
  
  delata_k_beta = rbinom(M,1,0.5); # [beta gamma rho]
  delata_k_gamma = rbinom(M,1,0.5);
  delata_k_rho = rbinom(M,1,0.5);
  
  delata_k = matrix(0,M,3)
  delata_k[,1] = delata_k_beta
  delata_k[,2] = delata_k_gamma
  delata_k[,3] = delata_k_rho
  delata_k[delata_k == 0] = -1
  
  c_delat_k = matrix(0,M,3)
  c_delat_k[,1] = c_m*delata_k[,1]
  c_delat_k[,2] = c_m*delata_k[,2]
  c_delat_k[,3] = c_m*delata_k[,3]
  
  # pre-allocate matricies 
  delta_G_k = matrix(0,3,1)
  c_delta_m = matrix(0,1,3)
  G_m_theta = matrix(0,3,1)
  
  for(m in 1:M){
    print(m)
    # Set J
    J = J_seq[m];
    
    theta_F = matrix(0,J,3)
    theta_P = matrix(0,J,3)
    
    # set marix for weigths
    
    w_m = matrix(0, N, J)
    
    # set indices
    index  = 1:J;
    
    # draw parameter perturbations
    Z_beta = rtruncnorm(J,-limit,limit, pd_beta_mu, pd_beta_sigma)
    Z_gamma = rtruncnorm(J,-limit,limit, pd_gamma_mu, pd_gamma_sigma) 
    Z_rho = rtruncnorm(J,-limit,limit, pd_rho_transformed_mu, pd_rho_transformed_sigma) 
    
    # Initilize parameters
    tau_Z_beta = tau_seq[m]*Z_beta
    tau_Z_gamma = tau_seq[m]*Z_gamma
    tau_Z_rho = tau_seq[m]*Z_rho
    
    theta_F[,1] = theta_vec[m,1] + tau_Z_beta
    theta_F[,2] = theta_vec[m,2] + tau_Z_gamma
    theta_F[,3] = theta_vec[m,3] + tau_Z_rho
    
    # Initilize states
    X_F = round( rmvnorm(J,c(30000, 800, 470000, 400), diag(c(50,30,500,10)) ) ) # X_F = [S I R H]
    
    popsize = rowSums(X_F[,1:3]); 
    
    # calc theta_bar_0
    theta_bar_0 = colMeans(theta_F);
    theta_bar_0_add = colMeans(theta_F) + c_delat_k[m,]
    theta_bar_0_sub = colMeans(theta_F) - c_delat_k[m,]
    
    for(n in 1:N){
      # Perturb parameters
      Z_beta = rtruncnorm(J,-limit,limit, pd_beta_mu, pd_beta_sigma)
      Z_gamma = rtruncnorm(J,-limit,limit, pd_gamma_mu, pd_gamma_sigma) 
      Z_rho = rtruncnorm(J,-limit,limit, pd_rho_transformed_mu, pd_rho_transformed_sigma) 
      
      theta_P[,1] = theta_F[,1] + sigma_seq[m]*Z_beta
      theta_P[,2] = theta_F[,2] + sigma_seq[m]*Z_gamma
      theta_P[,3] = theta_F[,3] + sigma_seq[m]*Z_rho
      
      
      # Simulate prediction particles
      X_P = simulatior_X(theta_P,X_F, popsize) # need to fix this
      
      # Evaluate weigths
      w = evaluator_Y_X(Y_hat[n], X_P, theta_P)
      prob_resample = w/sum(w) # normalize 
      # store wegiths
      
      w_m[n,] = w 
      
      # Systematic resampling of indecies 
      K = sample(index, J, replace = TRUE, prob = prob_resample)
      
      # Resaple 
      X_F = X_P[K,]
      theta_F = theta_P[K,]
      
      # calc mean
      theta_bar_vec[n, ] = colMeans(theta_F);
      theta_bar_vec_add[n,] = colMeans(theta_F) + c_delat_k[m,];
      theta_bar_vec_sub[n,] = colMeans(theta_F) - c_delat_k[m,];
      
      # calc variance 
      if(n == 1){
        theta_var_vec[n,] = (1/(J-1))*c(sum((theta_P[,1] - theta_bar_0[1])^2), 
                                        sum((theta_P[,2] - theta_bar_0[2])^2),
                                        sum((theta_P[,3] - theta_bar_0[3])^2));
        
        theta_var_vec_add[n,] = (1/(J-1))*c(sum((theta_P[,1] + c_delat_k[m,1] - theta_bar_0_add[1])^2),
                                             sum((theta_P[,2] + c_delat_k[m,2] - theta_bar_0_add[2])^2) ,
                                             sum((theta_P[,3] + c_delat_k[m,3] - theta_bar_0_add[3])^2));
        
        
        theta_var_vec_sub[n,] = (1/(J-1))*c(sum((theta_P[,1] - c_delat_k[m,1] - theta_bar_0_add[1])^2),
                                            sum((theta_P[,2] - c_delat_k[m,2] - theta_bar_0_add[2])^2),
                                            sum((theta_P[,3] - c_delat_k[m,3] - theta_bar_0_add[3])^2));
        
      }else{
        theta_var_vec[n,] = (1/(J-1))*c(sum((theta_P[,1] - theta_bar_vec[n-1,1])^2),  
                                        sum((theta_P[,2] - theta_bar_vec[n-1,2])^2), 
                                        sum((theta_P[,3] - theta_bar_vec[n-1,3])^2));
        
        theta_var_vec_add[n,] = (1/(J-1))*c(sum((theta_P[,1] + c_delat_k[m,1] - theta_bar_vec_add[n-1,1])^2),
                                             sum((theta_P[,2] + c_delat_k[m,2] - theta_bar_vec_add[n-1,2])^2),
                                             sum((theta_P[,3] + c_delat_k[m,3] - theta_bar_vec_add[n-1,3])^2));
        
        
        theta_var_vec_sub[n,] = (1/(J-1))*c(sum((theta_P[,1] - c_delat_k[m,1] - theta_bar_vec_add[n-1,1])^2),
                                            sum((theta_P[,2] - c_delat_k[m,2] - theta_bar_vec_add[n-1,2])^2),
                                            sum((theta_P[,3] - c_delat_k[m,3] - theta_bar_vec_add[n-1,3])^2));
        
      }
      
    }
    # calc log lik
    log_lik_vec[m] = sum(log(rowMeans(w_m)));
    
    # update 
    theta_var_vec = theta_var_vec^(-1);
    end_val = nrow(theta_bar_vec)
    G_m_theta[,1] = c(sum( theta_var_vec[,1]*(theta_bar_vec[1:N,1] - c(theta_bar_0[1], theta_bar_vec[1:(end_val-1),1]))),
                  sum( theta_var_vec[,2]*(theta_bar_vec[1:N,2] - c(theta_bar_0[2], theta_bar_vec[1:(end_val-1),2]))),
                  sum( theta_var_vec[,3]*(theta_bar_vec[1:N,3] - c(theta_bar_0[3], theta_bar_vec[1:(end_val-1),3]))))
    
    
    # calc delta G_k 
    theta_var_vec_add_beta = theta_var_vec_add[,1]^(-1);
    theta_var_vec_add_gamma = theta_var_vec_add[,2]^(-1);
    theta_var_vec_add_rho = theta_var_vec_add[,3]^(-1);
    
    theta_var_vec_sub_beta = theta_var_vec_sub[,1]^(-1);
    theta_var_vec_sub_gamma = theta_var_vec_sub[,2]^(-1);
    theta_var_vec_sub_rho = theta_var_vec_sub[,3]^(-1);
    
    end_val = nrow(theta_bar_vec)
    
    delta_G_k[,1] = c( sum( theta_var_vec_add_beta*(theta_bar_vec[1:N,1] - c( theta_bar_0_add[1], theta_bar_vec[1:(end_val-1),1] ) )) - 
                 sum( theta_var_vec_sub_beta*(theta_bar_vec[1:N,1] - c( theta_bar_0_sub[1], theta_bar_vec[1:(end_val-1),1]))), 
                 sum( theta_var_vec_add_gamma*(theta_bar_vec[1:N,2] - c(theta_bar_0_add[2], theta_bar_vec[1:(end_val-1),2] ) ) ) - 
                 sum( theta_var_vec_sub_gamma*(theta_bar_vec[1:N,2] - c(theta_bar_0_sub[2], theta_bar_vec[1:(end_val-1),2] ) ) ),
                 sum( theta_var_vec_add_rho*(theta_bar_vec[1:N,3] - c(theta_bar_0_add[3], theta_bar_vec[1:(end_val-1),3] ) ) ) - 
                 sum( theta_var_vec_sub_rho*(theta_bar_vec[1:N,3] - c(theta_bar_0_sub[3], theta_bar_vec[1:(end_val-1),3] ) ) ) );
    
    c_delta_m = c_delat_k[m,]
    
    H_hat_k = 0.5*(  delta_G_k%*%((2*c_delta_m)^(-1)) + t(  delta_G_k%*%((2*c_delta_m^(-1)) ) ) );
    
    
    if( m == 1 ) {
      H_bar_k = (1/(m+1))*H_hat_k;
    }else{
      H_bar_k = (m/(m+1))*H_bar_k + (1/(m+1))*H_hat_k;
    }
    
    #H_bar_bar_k = Re(sqrtm(H_bar_k*H_bar_k)) + matrix(0,3,3);
    
    d_k_vec = 10000*seq(1,M,1);
    
    d_k = d_k_vec[m];
    
    #H_bar_bar_k = real(sqrtm(H_bar_k*H_bar_k)) + d_k*eye(3);
    H_bar_bar_k = H_bar_k + d_k*diag(3);
    
    
    # update theta
    if( m >= burn_in) {
      theta_vec[m+1,] = theta_vec[m,] - a_m[m,]*t(ginv(H_bar_bar_k) %*% G_m_theta)
                                                     
    } else {
      theta_vec[m+1,] = theta_vec[m,];
    } 
    
  }
  
  
  # create output list
  output_list = list('theta_AIF' = theta_vec, 'loglik' = log_lik_vec)
  
  return(output_list)
}

# simulate_process
simulate_process <- function(N, nbr_of_steps){
  #library("truncnorm", lib.loc="~/R/win-library/3.2")
  theta = c(400, 26, 1/50, 0.1)
  X = c(30000, 800, 496200, 400)
  popsize = sum(X[1:3])
  Y = 50
  dt = 1/52/nbr_of_steps 
  
  Y_val <- rep(0, N)
  X_val <- matrix(0, N, 4)
  
  X_val[1,] = X
  Y_val[1] = Y
  
  
  for(n in 2:N) 
  {
    X = X_val[n-1,]
    for(i in 1:nbr_of_steps) 
    {
      
      if (i == 1) {
        P = sum(X[1:3])
        S = X[1]
        I = X[2]
        R = X[3]
        H = 0
      }else{
        P = sum(X[1:3])
        S = X[1]
        I = X[2]
        R = X[3]
        H = X[4]
      }
      
      beta = theta[1]
      gamma = theta[2]
      mu = theta[3]
      rho = theta[4]
      
      lambda_t = beta*I/popsize
      
      d_N0 = rpois(1, mu*popsize*dt)
      
      p_SI = (lambda_t/(lambda_t + mu))*(-expm1(-(lambda_t + mu)*dt))
      p_S = (mu/(lambda_t + mu))*(-expm1(-(lambda_t + mu)*dt))
      dNS = rmultinom(1,S,c(p_SI, p_S, 1-p_SI-p_S))
      dN_SI = dNS[1]
      dN_S = dNS[2]
      S_dN_SI_dN_S = dNS[3]
      S_temp = S + d_N0 - dN_SI - dN_S
      
      
      p_IR = (gamma/(gamma + mu))*(-expm1(-(gamma + mu)*dt))
      p_I = (mu/(gamma + mu))*(-expm1(-(gamma + mu)*dt))
      dNI = rmultinom(1,I,c(p_IR, p_I, 1-p_IR-p_I))
      dN_IR = dNI[1]
      dN_I = dNI[2]
      I_dN_IR_dN_I = dNI[3]
      I_temp = I + dN_SI - dN_IR - dN_I
      
      p_R = -expm1(-(mu)*dt)
      dNR = rmultinom(1 , R, c(p_R, 1- p_R))
      dN_R = dNR[1]
      R_dN_R = dNR[2]
      R_temp = R + dN_IR - dN_R
      H_temp = H + dN_IR
      
      
      X[1] = S_temp
      X[2] = I_temp
      X[3] = R_temp
      X[4] = H_temp
    }
    
    
    X_val[n, ] = X
    Y_val[n] = round(rho*X[4])
    if( Y_val[n] <= 0) 
    {
      Y_val[n] = 1
    }
    
  }
  data_list <- list("X_process" = X_val, "Y_process" = Y_val)
  return(data_list)
}

# simulatior_X
simulatior_X <- function(theta_P,X_F, popsize){
  
  X_P = matrix(0, nrow(X_F), ncol(X_F))
  nbr_of_steps = 5
  dt = 1/52/nbr_of_steps # migth need to use this somewhere...
  
  for( i in 1:nbr_of_steps ){
  
    if( i == 1 ){ 
      P = rowSums(X_F[,1:3]);
      S = X_F[,1]; # X_F = [S I R H]
      I = X_F[,2];
      R = X_F[,3];
      H = rep(0,nrow(X_F))
    }else{
      P = rowSums(X_P[,1:3]);
      S = X_P[,1];
      I = X_P[,2];
      R = X_P[,3];
      H = X_P[,4];
    }
    
    beta = theta_P[,1]; # theat_P = [beta gamma rho]
    gamma = theta_P[,2];
    mu = 1/50*rep(1, nrow(theta_P));
    
    lambta_t = beta*I/popsize;
    
    dNS = matrix(0, nrow(theta_P), 3)
    dNI = matrix(0, nrow(theta_P), 3)
    dNR = matrix(0, nrow(theta_P), 2)
    
    for( i in 1:nrow(theta_P) ) {
      
      lambta_t_temp = lambta_t[i]
      mu_temp = mu[i]
      gamma_temp = gamma[i]
      popsize_temp = popsize[i]
      
      d_N0 = rpois(1, mu_temp*popsize_temp*dt)
      
      
      p_SI = (lambta_t_temp/(lambta_t_temp + mu_temp))*(-expm1(-(lambta_t_temp + mu_temp)*dt))
      p_S = (mu_temp/(lambta_t_temp + mu_temp))*(-expm1(-(lambta_t_temp + mu_temp)*dt))
      
      if(p_SI < 0 | is.nan(p_SI) ){
        p_SI <- 10^-6
        print('Bad probability value')
      } 
      if(p_S < 0 | is.nan(p_S) ){
        p_S <- 10^-6
        print('Bad probability value')
      }  
      
      dNS = rmultinom(1,S[i],c(p_SI, p_S, 1-p_SI-p_S))
      dN_SI = dNS[1]
      dN_S = dNS[2]
      S_dN_SI_dN_S = dNS[3]
      S_temp = S[i] + d_N0 - dN_SI - dN_S
      
      
      p_IR = (gamma_temp/(gamma_temp + mu_temp))*(-expm1(-(gamma_temp + mu_temp)*dt))
      p_I = (mu_temp/(gamma_temp + mu_temp))*(-expm1(-(gamma_temp + mu_temp)*dt))
      
      if(p_IR < 0 | is.nan(p_IR) ){
        p_IR <- 10^-6
        print('Bad probability value')
      } 
      if(p_I < 0 | is.nan(p_I) ){
        p_I <- 10^-6
        print('Bad probability value')
      }  
      
      dNI = rmultinom(1,I[i],c(p_IR, p_I, 1-p_IR-p_I))
      dN_IR = dNI[1]
      dN_I = dNI[2]
      I_dN_IR_dN_I = dNI[3]
      I_temp = I[i] + dN_SI - dN_IR - dN_I
    
      
      p_R = -expm1(-(mu_temp)*dt)
      dNR = rmultinom(1 , R[i], c(p_R, 1 - p_R))
      dN_R = dNR[1]
      R_dN_R = dNR[2]
      R_temp = R[i] + dN_IR - dN_R
      H_temp = H[i] + dN_IR
      
      
      X_P[i,1] = S_temp
      X_P[i,2] = I_temp
      X_P[i,3] = R_temp
      X_P[i,4] = H_temp
      
    }
  
  }
  return(X_P)
}

# evaluator_Y_X
evaluator_Y_X <- function( Y_hat, X_P, theta_P ){   
  theta <- 100
  rho <- pnorm(theta_P[,3])
  H <- X_P[,4]
  mu = rho*H
  w = rep(0,nrow(theta_P))
  
  for(i in 1:nrow(theta_P)){
    w[i] = ( gamma( Y_hat + theta ) / ( gamma(theta)*factorial(Y_hat) ) )*( (  mu[i] /  ( mu[i] + theta ) )^(Y_hat) )*( ( 1 + ( mu[i]/theta ) )^(-theta) );
  }
  if( sum(is.nan(w)) > 0) {
    w[isnan(w)] <- 0;
    print('Weigths are NaN')
  }
  return(w)

}


