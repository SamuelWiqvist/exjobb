%% Set paramerets
i = 2
M = 100;
N = 500;
delta = 2;
K = 10;
m = 1:M;

[X,Y_hat] = generate_data(N);

% theta0,Y_hat, M, N,J_seq,tau_seq,a_seq,sigma_seq,Sigma

J = 1000;
theta0 = -0.45; %-1 + (1+1)*rand(1,1); % -0.6407 seems to work well
J_seq_fast = 50.*m.^(2);
J_seq = 50.*m;%.*m; % m.^(delta + 0.5); 

J0 = 1000;
J1 = 200;
J2 = 100;
J_seq_d_g = J0*exp(-(m-1)) + J1 + J2*m.^(0.6);

J_seq_c = 1000*ones(M,1);
tau_seq = sqrt(m.^(-1));
sigma_seq = sqrt(m.^(-delta - 0.5));

alpha = 1;
a = 0.1;
A = 0.05*M;
a_seq = (m + 1 + A).^(-alpha).*a;
Sigma = 0.2;

factor = 0.2;
sigma_m_2_1 = (1:M).^(-0.6)*factor;
sigma_m_2_2 = (1:M).^(-0.95)*factor;
sigma_a1 = 0.99.^(0:M-1)*factor;

theta_all_d_g_J = zeros(M,K);
theta_all_changefastJ = zeros(M,K);
theta_all_constantJ = zeros(M,K);
theta_all_constantJ_opt = zeros(M,K);

theta_all_IF1 = zeros(M,K);
theta_all_hessia = zeros(M,K);
matrix_conv = zeros(M,K);
theta_all_kestans = zeros(M,K);

loglik_all_IF1 = zeros(M,K);
loglik_all_hessian = zeros(M,K);
loglik_all_kestans = zeros(M,K);

theta_all_IF2 = zeros(M,K);

% parameters for AIF
burn_in = 10;
M_A = M + burn_in;
m_A = 1:M_A;

a_m = 50./((m_A + 1 + 0.05*M).^(0.9));
c_m = (std(Y_hat)*2)./((m_A + 1).^(0.2));
J_seq_A = 50*ones(M_A,1);  %.*m; % m.^(delta + 0.5); 
tau_seq_A = sqrt(m_A.^(-1));
sigma_seq_A = sqrt(m_A.^(-delta - 0.5));
a_k_burn_in = a_seq(1:burn_in);

if i == 1
    a_k = 1./((m + 0.1*M).^(1)); % asymtotcally optimal parameters are realy bad... 
    c_k = (std(Y_hat)*2)./((m + 1).^(1/6));
    a_k_v2 = 1./((m_A + 0.1*M_A).^(1)); % asymtotcally optimal parameters are realy bad... 
    c_k_v2 = (std(Y_hat)*2)./((m_A + 1).^(1/6));

elseif i == 2
    a_k = 1./((m + 0.1*M).^(0.7)); % asymtotcally optimal parameters are realy bad... 
    c_k = (std(Y_hat)*5)./((m + 1).^(0.3));
    a_k_v2 = 1./((m_A + 0.1*M).^(0.7)); % asymtotcally optimal parameters are realy bad... 
    c_k_v2 = (std(Y_hat)*2)./((m_A + 1).^(0.3));

else
    a_k = 1./((m + 0.1*M).^(0.602)); % asymtotcally optimal parameters are realy bad... 
    c_k = (std(Y_hat)*2)./((m + 1).^(0.101));
    a_k_v2 = 1./((m_A + 0.1*M).^(0.602)); % asymtotcally optimal parameters are realy bad... 
    c_k_v2 = (std(Y_hat)*2)./((m_A + 1).^(0.101));

end


% calc theta_ML
theta_all_ML = ML_est(Y_hat,theta0)

%%
% Run algorithms 
for k = 1:K
    k
    
    nbr_of_restarts_IF1 = 0;
    while(true)
        try
            [theta_obs_IF1, loglikIF1] = IF1_v2(theta0,Y_hat, M, N,J_seq,tau_seq,a_seq,sigma_seq,Sigma);
            if exist('theta_obs_IF1', 'var')
                'Done: IF1_v2'
                break
            end
        catch ME
            switch ME.identifier
                case 'stats:datasample:InvalidWeights'
                    nbr_of_restarts_IF1 = nbr_of_restarts_IF1+1;
                    warning('Restart');
            end
        end
    end
    
    nbr_of_restarts_Hessia = 0;
    while(true)
        try
            [theta_obs_hessia, loglik_hessia, H_vec, nbr_accpted]= IF1_A_Hessian_approx(theta0,Y_hat, M, N,J_seq,tau_seq,sigma_seq,Sigma, a_k,c_k);
            if exist('theta_obs_hessia', 'var')
                'Done: IF1_A_Hessian_approx'
                break
            end
        catch ME
            switch ME.identifier
                case 'stats:datasample:InvalidWeights'
                    nbr_of_restarts_Hessia = nbr_of_restarts_Hessia+1;
                    warning('Restart');
            end
        end
    end
    

    
%     nbr_of_restarts_kestans = 0;
%     while(true)
%         try
%             [theta_obs_kestens, loglikkestens,k_hat] = IF1_A_kestens(theta0,Y_hat, M, N,J_seq,tau_seq,a_seq,sigma_seq,Sigma);
%             if exist('theta_obs_kestens', 'var')
%                 'Done: IF1_A_kestens'
%                 break
%             end
%         catch ME
%             switch ME.identifier
%                 case 'stats:datasample:InvalidWeights'
%                     nbr_of_restarts_kestans = nbr_of_restarts_kestans+1;
%                     warning('Restart');
%             end
%         end
%     end
    
    
    theta0_m = -1 + (1+1)*rand(J,1);
    theta_obs_IF2 = IF2(theta0_m,Y_hat, M, sigma_a1,J,N);
    theta_all_IF2(:,k) = mean(theta_obs_IF2(1:end-1,:),2);
    
    theta_all_hessia(:,k) = theta_obs_hessia;
    theta_all_IF1(:,k) = theta_obs_IF1;
    %theta_all_kestans(:,k) = theta_obs_kestens;
    
    loglik_all_IF1(:,k) = loglikIF1;
    loglik_all_hessian(:,k) = loglik_hessia;
    %loglik_all_kestans(:,k) = loglikkestens;
    
    matrix_conv(:,k) =  H_vec;

    
end


%%
theta_ML = mean(theta_all_ML);

% RMSE
RMSE_IF1 = sqrt(mean((theta_all_IF1 - theta_ML).^2,2));
RMSE_hessian = sqrt(mean((theta_all_hessia - theta_ML).^2,2));
RMSE_IF2 = sqrt(mean((theta_all_IF2 - theta_ML).^2,2));

figure
loglog(10*m.^(-2), '--')
hold on
loglog(10*m.^(-1), '--')
loglog(10*sqrt(m).^(-1), '--')
loglog(RMSE_IF1)
loglog(RMSE_hessian)
loglog(RMSE_IF2)
xlabel('Iteration')
ylabel('RMSE')
hleg1 = legend('O(1/m^2)', 'O(1/m)', 'O(1/sqrt(m))', 'IF1', 'AIF', 'IF2', 'Location','southwest');
set(hleg1,'FontSize',10)


%% Plot loglik
conv_phase_end = 300;
d = size(loglik_all_IF1);
d_col = d(2);

figure

subplot(211)
plot(loglik_all_IF1(1:conv_phase_end,:))
title('IF1')
ylabel('loglik')

subplot(212)
plot(loglik_all_hessian(1:conv_phase_end,:))
title('AIF1')
ylabel('loglik')


xlabel('Iteration')
ylabel('loglik')

figure

subplot(211)
plot([zeros(300,d_col); loglik_all_IF1(conv_phase_end+1:end,:)])
axis([300, 1000, -950, -900])
title('IF1')
ylabel('loglik')

subplot(212)
plot([zeros(300,d_col); loglik_all_hessian(conv_phase_end+1:end,:)]) 
axis([300, 1000, -950, -900])
title('AIF1')
ylabel('loglik')

xlabel('Iteration')
ylabel('loglik')

%% more plotting 
conv_phase_end = 300;
d = size(theta_all_IF1);
d_col = d(2);

theta_ture = 0.8*ones(length(theta_all_IF1),1);


figure

subplot(211)
plot(theta_all_IF1(1:conv_phase_end,:))
hold on
plot(theta_ture(1:conv_phase_end), 'k--')
title('IF1')
ylabel('\theta')

subplot(212)
plot(theta_all_hessia(1:conv_phase_end,:))
hold on
plot(theta_ture(1:conv_phase_end), 'k--')
title('AIF1')
ylabel('\theta')
xlabel('Iteration')


figure

subplot(211)
plot([zeros(300,d_col); theta_all_IF1(conv_phase_end+1:end,:)])
hold on
plot([zeros(300,1); theta_ture(conv_phase_end+1:end)], 'k--')
axis([300, 1000, 0.78 1])
title('IF1')
ylabel('\theta')

subplot(212)
plot([zeros(300,d_col); theta_all_hessia(conv_phase_end+1:end,:)])
hold on
plot([zeros(300,1); theta_ture(conv_phase_end+1:end)], 'k--')
axis([300, 1000, 0.78 1])
title('AIF1')
ylabel('\theta')
xlabel('Iteration')


figure
plot(matrix_conv, '--')
ylabel('Hessian ')
xlabel('Iteration')
