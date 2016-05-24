# General parameters
rm(list = ls())
source('~/R/functions.R')

M = 50;
N = 500;
delta = 2;
K = 3;
m = 1:M;



# Generate data
ptm <- proc.time()
nbr_of_steps = 20
data <- simulate_process(N,nbr_of_steps)
runtime <- proc.time() - ptm

Y_hat = data$Y_process

# Plot Data
time = seq(0, 1/52*N, 1/52);
time = time[1:N];



par(mfrow=c(2,2))
plot(time,data$X_process[,1], main="S")
plot(time,data$X_process[,2], main="I")
plot(time,data$X_process[,3], main="R")
plot(time,data$X_process[,4], main="H")

par(mfrow=c(1,1))
plot(time,data$Y_process, main="Cases")

write.table(data$Y_process, "./cases.txt", sep = " ", col.names = FALSE, row.names = FALSE)
write.table(data$X_process, "./X_process.txt", sep = " ", col.names = FALSE, row.names = FALSE)

# start values 
beta_0 = (500 - 300)*runif(1) + 300
gamma_0 = (40-10)*runif(1) + 10
rho_0 = (0.5-0)*runif(1) + 0
theta_0 = c(beta_0, gamma_0,  qnorm(rho_0))  # [beta gamma  rho]


# set parameters 
J_seq = 50*(1:M)#50*(1:M) #50.*m;%.*m; % m.^(delta + 0.5);  50*ones(M,1);%

tau_seq = sqrt(m^(-1))
sigma_seq = sqrt(m^(-delta - 0.5))

alpha = 1
a = 0.02
A = 0.05*M
a_seq = matrix(0,M,3)
a_seq[,1] = ((m + 1 + A)^(-alpha))*20
a_seq[,2] = ((m + 1 + A)^(-alpha))*1
a_seq[,3] = ((m + 1 + A)^(-alpha))*0.01

##
Sigma = diag(3)
Sigma[1,1] = 20
Sigma[2,2] = 0.2
Sigma[3,3] = 0.2

# parameters IF2

# start values IF2  

J = 500
theta_0_m = matrix(0,J,3) # [beta gamma  rho]
theta_0_m[,1] = (500 - 300)*runif(J) + 300
theta_0_m[,2] = (40-10)*runif(J) + 10
theta_0_m[,3] = (0.5-0)*runif(J) + 0

factor = 1
sigma_m = 0.99^(0:M-1)*factor

#parameters for AIF
burn_in = 1;
M_A = M + burn_in;
m_A = 1:M_A;

a_m = matrix(0,M_A,3)
a_m[,1] = 10/((m_A + 1 + 0.05*M)^(0.9))
a_m[,2] = 5/((m_A + 1 + 0.05*M)^(0.9))
a_m[,3] = 0.5/((m_A + 1 + 0.05*M)^(0.9))


c_m = (sd(Y_hat)*2)/((m_A + 1)^(0.2))

J_seq_A = 50*(1:M_A);  #.*m; % m.^(delta + 0.5); 
tau_seq_A = sqrt(m_A^(-1))
sigma_seq_A = sqrt(m_A^(-delta - 0.5))

# est par with IF1

ptm <- proc.time()
par_est_IF1 <- IF1_SIR(theta_0,Y_hat, M, 500,J_seq,tau_seq,a_seq,sigma_seq,Sigma)
runtime <- proc.time() - ptm

theta_obs = par_est_IF1$theta_IF1

par(mfrow=c(3,1))
plot(theta_obs[,1], main="beta") # [beta gamma rho]
plot(theta_obs[,2], main="gamma")
plot(pnorm(theta_obs[,3]), main="rho")

# est par with IF2
par_est_IF2 <- IF2_SIR(theta_0_m,Y_hat, M, sigma_m,J,N)

beta_IF2 = rowMeans(par_est_IF2[,1:J]);
gamma_IF2 = rowMeans(par_est_IF2[,J+1:2*J]);
rho_IF2 = pnorm(rowMeans((par_est_IF2[,(2*J+1):(3*J)])));

par(mfrow=c(3,1))
plot(beta_IF2, main="beta") # [beta gamma rho]
plot(gamma_IF2, main="gamma")
plot(rho_IF2, main="rho")

# est par with AIF

ptm <- proc.time()
par_est_AIF <- AIF_SIR(theta_0,Y_hat, M_A, 500,J_seq_A,tau_seq_A,sigma_seq_A,Sigma,a_m,c_m,burn_in)
runtime <- proc.time() - ptm


theta_obs_AIF = par_est_AIF$theta_AIF

par(mfrow=c(3,1))
plot(theta_obs_AIF[,1], main="beta") # [beta gamma rho]
plot(theta_obs_AIF[,2], main="gamma")
plot(pnorm(theta_obs_AIF[,3]), main="rho")

# Run all algorithms 

theta_IF1 = matrix(0,M+1,K*3);
theta_IF2 = matrix(0,M,K*3);
theta_AIF = matrix(0,M_A+1,K*3);

stor_m_IF2 = matrix(0,M,3)
start_ind = seq(1,K*3,3)

for (i in 1:K){
  par_est_IF1 <- IF1_SIR(theta_0,Y_hat, M, 500,J_seq,tau_seq,a_seq,sigma_seq,Sigma)
  
  par_est_IF2 <- IF2_SIR(theta_0_m,Y_hat, M, sigma_m,J,N)
  
  beta_IF2 = rowMeans(par_est_IF2[,1:J]);
  gamma_IF2 = rowMeans(par_est_IF2[,J+1:2*J]);
  rho_IF2 = rowMeans((par_est_IF2[,(2*J+1):(3*J)]));
  
  stor_m_IF2[,1] = beta_IF2
  stor_m_IF2[,2] = gamma_IF2
  stor_m_IF2[,3] = rho_IF2
  
  par_est_AIF <- AIF_SIR(theta_0,Y_hat, M_A, 500,J_seq_A,tau_seq_A,sigma_seq_A,Sigma,a_m,c_m,burn_in)

  theta_IF1[,start_ind[i]:(start_ind[i]+2)] = par_est_IF1$theta_IF1

  theta_IF2[,start_ind[i]:(start_ind[i]+2)] = stor_m_IF2
  theta_AIF[,start_ind[i]:(start_ind[i]+2)] = par_est_AIF$theta_AIF
}

theta_IF1 = theta_IF1[1:(M+1-1),]
theta_AIF = theta_AIF[(burn_in+1):(M_A+1),]

write.table(theta_IF1, "./data_IF1.txt", sep = " ", col.names = FALSE, row.names = FALSE)
write.table(theta_IF2, "./data_IF2.txt", sep = " ", col.names = FALSE, row.names = FALSE)
write.table(theta_AIF, "./data_AIF.txt", sep = " ", col.names = FALSE, row.names = FALSE)

