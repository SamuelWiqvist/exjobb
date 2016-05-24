%% load data
data = dlmread('data_g.txt');
Y_hat = data(:,2);
plot(Y_hat)
%% simulate data 
[X,Y_hat] = generate_data(1,theta_true, 1, 100);

%%
N = 100;
M = 200;
burn_in = 1;
M_A = M + burn_in;
m = 1:M;
m_A = 1:M_A;
K = 5;
k_true = 1;
delta = 2;
theta_true = [0.1,1,0.1,0.1]; % [r k sigma tau] 
Sigma = 0.02*eye(3);

r_0 = rand(1) ;
sigma_0 = rand(1);
tau_0 = rand(1);
theta0 = [r_0 sigma_0 tau_0];

J_seq = 50.*m;%.*m; % m.^(delta + 0.5); 
 
tau_seq = sqrt(m.^(-1));
sigma_seq = sqrt(m.^(-delta - 0.5));

alpha = 1;
a = 0.02;
A = 0.05*M;
a_seq = (m + 1 + A).^(-alpha).*a;
plot(exp(Y_hat))

theta_est = zeros(M,K*3);
loglik_est = zeros(M,K);
theta_A = zeros(M_A,K*3);
loglik_A = zeros(M_A,K);
theta_IF2 = zeros(M,K*3);


theta_est_var = zeros(M,K*2);
loglik_est_var = zeros(M,K);
theta_A_var = zeros(M_A,K*2);
loglik_A_var = zeros(M_A,K);
theta_IF2_var = zeros(M,K*2);

% parameter for the adaptive algorithm
a_m = 50./((m_A + 1 + 0.05*M).^(0.9));
c_m = (std(exp(Y_hat))*2)./((m_A + 1).^(0.2));
J_seq_A = 50.*m_A;  %.*m; % m.^(delta + 0.5); 
tau_seq_A = sqrt(m_A.^(-1));
sigma_seq_A = sqrt(m_A.^(-delta - 0.5));

% parameters for the IF2 algorithm
J = 1000;
r_0_IF2 = rand(1000,1);
sigma_0_IF2 = rand(1000,1);
tau_0_IF2 = rand(1000,1);
theta0_IF2 = [r_0_IF2, sigma_0_IF2, tau_0_IF2];
theta0_IF2_est_var = [sigma_0_IF2, tau_0_IF2];
factor = 0.001;
sigma_a1 = 0.99.^(0:(M-1))*factor;
plot(sigma_a1)
%% Est r sigma tau with AIF 
[theta_obs_A, log_lik_vec, nbr_not_a, eig_vec, norm_vec] = IF1_gompertz_v2_Hessian(log(theta0),exp(Y_hat), M_A, N,J_seq_A,tau_seq_A,sigma_seq_A,Sigma,k_true, true, a_m,c_m,burn_in);

exp(theta_obs_A(end,:))
figure
subplot(131)
plot(exp(theta_obs_A(burn_in:end-1,1)))
title('r')
subplot(132)
plot(exp(theta_obs_A(burn_in:end-1,2)))
title('\sigma')
subplot(133)
plot(exp(theta_obs_A(burn_in:end-1,3)))
title('\tau')
%% Est r sigma tau with IF1
[theta_obs_changeJ, loglikchangeJ] = IF1_gompertz_v2(log(theta0),exp(Y_hat), M, N,J_seq,tau_seq,a_seq,sigma_seq,Sigma,k_true, true);
exp(theta_obs_changeJ(end,:))
figure
subplot(131)
plot(exp(theta_obs_changeJ(:,1)))
subplot(132)
plot(exp(theta_obs_changeJ(:,2)))
title('IF1')
subplot(133)
plot(exp(theta_obs_changeJ(:,3)))

%% Est r sigma tau with IF2 
theta_obs_IF2 = IF2_gompertz(log(theta0_IF2),exp(Y_hat), M, sigma_a1,J,N,k_true);

r_IF2 = exp(mean(theta_obs_IF2(:,1:J),2));
sigma_IF2 = exp(mean(theta_obs_IF2(:,J+1:2*J),2));
tau_IF2 = exp(mean(theta_obs_IF2(:,2*J:end),2));

[r_IF2(end) sigma_IF2(end) tau_IF2(end)]

figure
subplot(131)
plot(r_IF2)
subplot(132)
plot(sigma_IF2)
subplot(133)
plot(tau_IF2)
