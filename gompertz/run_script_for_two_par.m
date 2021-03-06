%% load data
data = dlmread('data_g.txt');
Y_hat = data(:,2);
plot(Y_hat)

%%
N = 100;
M = 1000;
burn_in = 5;
M_A = M + burn_in;
m = 1:M;
m_A = 1:M_A;
K = 10;
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


% parameter for the adaptive algorithm
a_m = 50./((m_A + 1 + 0.05*M).^(0.9));
c_m = (std(exp(Y_hat))*10)./((m_A + 1).^(0.2));
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
factor = 0.0001;
sigma_a1 = 0.996.^(0:(M-1))*factor;
%sigma_a1 = m.^(-0.6)*factor
plot(sigma_a1)
%% Test with K est sigma,tau 

theta_est_var = zeros(M,K*2);
loglik_est_var = zeros(M,K);
theta_A_var = zeros(M_A,K*2);
loglik_A_var = zeros(M_A,K);
theta_IF2_var = zeros(M,K*2);

start_ind = 1:2:K*2;
for k = 1:K
    k
    nbr_of_restarts_changeJ = 0;
    while(true)
        try
            [theta_obs_changeJ, loglikchangeJ] = IF1_gompertz_v2_est_var(log(theta0(2:3)),exp(Y_hat), M, N,J_seq,tau_seq,a_seq,sigma_seq,Sigma,k_true, false);
            if exist('theta_obs_changeJ', 'var')
                %'Done: AIF'
                break
            end
        catch ME
            switch ME.identifier
                case 'stats:datasample:InvalidWeights'
                    nbr_of_restarts_changeJ = nbr_of_restarts_changeJ+1;
                    warning('Restart');
            end
        end
    end
    
    nbr_of_restarts_changeJ = 0;
    while(true)
        try
            [theta_obs_A, log_lik_vec, nbr_not_a, eig_vec, norm_vec] = IF1_gompertz_v2_Hessian_est_var(log(theta0(2:3)),exp(Y_hat), M_A, N,J_seq_A,tau_seq_A,sigma_seq_A,Sigma,k_true, false, a_m,c_m,burn_in);
            if exist('theta_obs_A', 'var')
                %'Done: AIF'
                break
            end
        catch ME
            switch ME.identifier
                case 'stats:datasample:InvalidWeights'
                    nbr_of_restarts_changeJ = nbr_of_restarts_changeJ+1;
                    warning('Restart');
            end
        end
    end
    
    theta_obs_IF2 = IF2_gompertz_est_var(log(theta0_IF2_est_var),exp(Y_hat), M, sigma_a1,J,N,k_true);


    sigma_IF2 = mean(theta_obs_IF2(:,1:J),2);
    tau_IF2 = mean(theta_obs_IF2(:,J+1:2*J),2);


    theta_est_var(:,start_ind(k):start_ind(k)+1) = theta_obs_changeJ;
    loglik_est_var(:,k) = loglikchangeJ;
    loglik_A_var(:,k) = log_lik_vec;
    theta_A_var(:,start_ind(k):start_ind(k)+1) = theta_obs_A;
    
    theta_IF2_var(:,start_ind(k):start_ind(k)+1) = [ sigma_IF2 tau_IF2];
end 

theta_est_var = exp(theta_est_var);
theta_A_var = exp(theta_A_var(burn_in+1:end,:));
theta_IF2_var = exp(theta_IF2_var);




%% plot RMSE est r,sigma,tau 

RMSE_IF_sigma = sqrt(mean((theta_est_var(:,start_ind) - 0.0322).^2,2));
RMSE_IF_tau = sqrt(mean((theta_est_var(:,start_ind+1) - 0.0694).^2,2));

RMSE_IF_sigma_A = sqrt(mean((theta_A_var(:,start_ind) - 0.0322).^2,2));
RMSE_IF_tau_A = sqrt(mean((theta_A_var(:,start_ind+1) - 0.0694).^2,2));

RMSE_IF_sigma_IF2 = sqrt(mean((theta_IF2_var(:,start_ind) - 0.0322).^2,2));
RMSE_IF_tau_IF2 = sqrt(mean((theta_IF2_var(:,start_ind+1) - 0.0694).^2,2));



figure
loglog(10*m.^(-1), '--')
hold on
loglog(10*sqrt(m).^(-1), '--')
loglog(RMSE_IF_sigma)
loglog(RMSE_IF_sigma_A)
loglog(RMSE_IF_sigma_IF2)
xlabel('Iteration')
ylabel('RMSE')
hleg1 = legend('O(1/m)', 'O(1/sqrt(m))', 'IF1', 'AIF1', 'IF2', 'Location','southwest')
set(hleg1,'FontSize',12)

figure
loglog(10*m.^(-1), '--')
hold on
loglog(10*sqrt(m).^(-1), '--')
loglog(RMSE_IF_tau)
loglog(RMSE_IF_tau_A)
loglog(RMSE_IF_tau_IF2)
xlabel('Iteration')
ylabel('RMSE')
hleg1 = legend('O(1/m)', 'O(1/sqrt(m))', 'IF1', 'AIF1', 'IF2', 'Location','southwest')
set(hleg1,'FontSize',12)




%% Plot convergence plots est r,sigma,tau 
theta_ture = 0.1*ones(length(theta_est_var),1);

% IF1
figure
subplot(121)
plot(theta_est_var(:,start_ind))
hold on
plot(theta_ture, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\sigma')
subplot(122)
plot(theta_est_var(:,start_ind+1))
hold on
plot(theta_ture, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\tau')

est_IF1 = [mean(theta_est_var(end,start_ind)), mean(theta_est_var(end,start_ind+1)) ] 

% IF2
figure
subplot(121)
plot(theta_IF2_var(:,start_ind))
hold on
plot(theta_ture, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\sigma')
subplot(122)
plot(theta_IF2_var(:,start_ind+1))
hold on
plot(theta_ture, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\tau')

est_IF2 = [mean(theta_IF2_var(end,start_ind)), mean(theta_IF2_var(end,start_ind+1)) ] 

theta_ture_a = 0.1*ones(1,length(theta_A_var));

% AIF
figure
subplot(121)
plot(theta_A_var(:,start_ind))
hold on
plot(theta_ture_a, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\sigma')
subplot(122)
plot(theta_A_var(:,start_ind+1))
hold on
plot(theta_ture_a, 'k--')
xlabel('Iteration')
ylabel('Value')
title('\tau')

est_AIF = [mean(theta_A_var(end,start_ind)), mean(theta_A_var(end,start_ind+1)) ] 



