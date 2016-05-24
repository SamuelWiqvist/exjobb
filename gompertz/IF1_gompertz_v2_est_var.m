function [theta_obs, log_lik_vec] = IF1_gompertz_v2_est_var(theta0_log,Y_hat, M, N,J_seq,tau_seq,a_seq,sigma_seq,Sigma,par_K, print_m)
    % Description
    % The IF2 algorithm based upon the algorithm described in
    % Ionides (2011)
    %
    % Input:
    % theta0 = start value 
    % Y_hat = observatrions 
    % M = nbr of iterations 
    % N = nbr of timestep 
    % J_seq = nbr of particles
    % tau_seq = cooling factor for theta
    % sigma_seq = cooling factor for theta
    % a_seq = cooling factor for the updating formula
    % Sigma = covariance matrix for the distribtuion k 
    %
    % Output:
    % theta_obs = the estimated values of theta for each M-iteration
    % log_lik_vec = estiamtion of the log-lik
    
    % Print msg 
    disp('Starting IF1_gompertz_v2')
    
    r = 0.1;
    
    log_lik_vec = zeros(M,1);
    
    theta_vec = zeros(M+1,2);
    theta_vec(1,:) = theta0_log;
    
    % set distribution for parameter perturbation     
    pd_sigma = makedist('Normal');
    limit = 1000;
    pd_sigma.mu = 0;
    pd_sigma.sigma = Sigma(2,2);
    kappa_sigma = truncate(pd_sigma,-limit,limit);
    
    pd_tau = makedist('Normal');
    limit = 1000;
    pd_tau.mu = 0;
    pd_tau.sigma = Sigma(3,3);
    kappa_tau = truncate(pd_tau,-limit,limit);
    
    theta_bar_vec_log = zeros(N,2);
    theta_var_vec_log = zeros(N,2);
    

    for m = 1:M
        if print_m == true
            m
        end
        % Set J
        J = J_seq(m);
        
        % set marix for weigths
        
        w_m = zeros(N,J);
        % set indices
        index  = 1:J;
        
        % draw parameter perturbations
        Z_sigma = random(kappa_sigma,J,1);
        Z_tau = random(kappa_tau,J,1);

        
        % Initilize parameters
        tau_Z = tau_seq(m).*[Z_sigma, Z_tau];
        theta_F_log = [theta_vec(m,1) + tau_Z(:,1), theta_vec(m,2) + tau_Z(:,2)];
        
        % Initilize states 
        X_F = rand(J,1);
        
        % calc theta_bar_0
        theta_bar_0_log = mean(theta_F_log);
        for n = 1:N
            
            % Perturb parameters
            Z_sigma = random(kappa_sigma,J,1);
            Z_tau = random(kappa_tau,J,1);
            
            theta_P_log = theta_F_log + sigma_seq(m).*[Z_sigma, Z_tau];
            
            % Simulate prdiction particles
            eps_temp = exp(normrnd(0,real(exp(theta_P_log(:,1)))));
            X_P = par_K.^(1-exp(-r)).*X_F.*(exp(-r)).*eps_temp;%normrnd(theta_P.*X_F,1);
            
            % Evaluate weigths
            w = normpdf(log(Y_hat(n)), log(X_P) , real(exp(theta_P_log(:,2)))); 
            prob_resample = w./sum(w); % normalize , set resampling probability 
            
            % store wegiths
            
            w_m(n,:) = w; 
            % Systematic resampling of indecies 
            k = datasample(index,J,'Replace',true, 'Weights', prob_resample);
            
            % Resaple 
            X_F = X_P(k);
            theta_F_log = theta_P_log(k,:);
            
            % calc mean
            theta_bar_vec_log(n,:) = mean(theta_F_log);
            
            theta_P = theta_P_log;
            if n == 1
                theta_var_vec_log(n,:) = [(1/(J-1))*sum((theta_P(:,1) - theta_bar_0_log(1)).^2), ...
                                          (1/(J-1))*sum((theta_P(:,2) - theta_bar_0_log(2)).^2)];
            else
                theta_var_vec_log(n,:) = [(1/(J-1))*sum((theta_P(:,1) - theta_bar_vec_log(n-1,1)).^2), ...
                                          (1/(J-1))*sum((theta_P(:,2) - theta_bar_vec_log(n-1,2)).^2)];
            end
        end
        
        % calc log lik
        log_lik_vec(m) = sum(log(mean(w_m,2)));
        
        % update 
        theta_var_vec_log = theta_var_vec_log.^(-1);
        
        G_m_theta = [sum( theta_var_vec_log(:,1).*(theta_bar_vec_log(1:N,1) - [theta_bar_0_log(1); theta_bar_vec_log(1:end-1,1)])), ...
                     sum( theta_var_vec_log(:,2).*(theta_bar_vec_log(1:N,2) - [theta_bar_0_log(2); theta_bar_vec_log(1:end-1,2)]))];
                                               
        theta_vec(m+1,:) = theta_vec(m,:) + a_seq(m)*G_m_theta;
        
        % accept update
        if abs(theta_vec(m+1,1) - theta_vec(m,1)) > 1 ||  abs(theta_vec(m+1,2) - theta_vec(m,2)) > 1 
            theta_vec(m+1,:) = theta_vec(m,:);
            print('Bad step')
        end
        
    end
   
    theta_obs =   theta_vec(1:end-1,:);    
end
