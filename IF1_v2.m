function [theta_obs, log_lik_vec] = IF1_v2(theta0,Y_hat, M, N,J_seq,tau_seq,a_seq,sigma_seq,Sigma)
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
    disp('Starting IF1_v2')
    
    log_lik_vec = zeros(M,1);
    
    theta_vec = zeros(M+1,1);
    theta_vec(1) = theta0;
    
    % set distribution for parameter perturbation 
    pd = makedist('Normal');
    pd.mu = 0;
    pd.sigma = Sigma;
    kappa = truncate(pd,-100,100);
    
    theta_bar_vec = zeros(N,1);
    theta_var_vec = zeros(N,1);
    
    
    for m = 1:M
        
        % Set J
        J = J_seq(m);
        
        % set marix for weigths
        
        w_m = zeros(N,J);
        % set indices
        index  = 1:J;
        
        % draw parameter perturbations
        Z = random(kappa,J,1);
        
        % Initilize parameters
        theta_F = theta_vec(m) + tau_seq(m).*Z;
        
        % Initilize states 
        X_F = normrnd(0, 1, J,1);
        
        % calc theta_bar_0
        theta_bar_0 = mean(theta_F);
        for n = 1:N
            % Perturb parameters
            Z = random(kappa,J,1);
            theta_P = theta_F + sigma_seq(m).*Z;
            
            % Simulate prdiction particles
            X_P = normrnd(theta_P.*X_F,1);
            
            % Evaluate weigths
            w = normpdf(Y_hat(n), X_P , 1); 
            prob_resampling = w./sum(w); % normalize 
            
            % store wegiths
            
            w_m(n,:) = w; 
            % Systematic resampling of indecies 
            k = datasample(index,J,'Replace',true, 'Weights', prob_resampling);
            
            % Resaple 
            X_F = X_P(k);
            theta_F = theta_P(k);
            
            % calc mean
            theta_bar_vec(n) = mean(theta_F);
            if n == 1
                theta_var_vec(n) = (1/(J-1))*sum((theta_P - theta_bar_0).^2);
            else
                theta_var_vec(n) = (1/(J-1))*sum((theta_P - theta_bar_vec(n-1)).^2);
            end
        end
        
        % calc log lik
        log_lik_vec(m) = sum(log(mean(w_m,2)));
        
        % update 
        theta_vec(m+1) = theta_vec(m) + a_seq(m)*sum( theta_var_vec.^(-1).*(theta_bar_vec(1:N) - [theta_bar_0; theta_bar_vec(1:end-1)]));
        
    end
   
    theta_obs =   theta_vec(1:end-1);    
end
