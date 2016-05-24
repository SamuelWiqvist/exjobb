function [theta_obs, log_lik_vec,H_vec, nbr_not_a] = IF1_A_Hessian_approx(theta0,Y_hat, M, N,J_seq,tau_seq,sigma_seq,Sigma, a_m,c_m)
    % Description: 
    % The IF1 algorithm with an adaptive updateing formula 
    % based upon the algorithm described in Ionides (2011) 
    % and (Spall 2000)
    %
    % Input:  
    % theta0 = start value 
    % Y_hat = observatrions 
    % M = nbr of iterations 
    % N = nbr of timestep 
    % J_seq = nbr of particles
    % tau_seq = cooling factor for theta
    % sigma_seq = cooling factor for theta
    % Sigma = covariance matrix for the distribtuion k 
    % a_m = scalar gain cooling coefficient 
    % c_m = parameter for pre-iterated estiamte of the Hessian 
    %
    % Output:
    % theta_obs = the estimated values of theta for each M-iteration
    % log_lik_vec = estiamtion of the log-lik
    % nbr_not_a = nbr of steps not accepted 
    
    % Print msg 
    disp('Starting IF1_A_Hessian_approx')
    
    % hessians
    H_vec = zeros(M,1);
    % set start value for the nbr of steps that are not accepted
    nbr_not_a = 1;
    
    log_lik_vec = zeros(M,1);
    
    theta_vec = zeros(M+1,1);
    theta_vec(1) = theta0;
    
    % set distribution for parameter perturbation 
    pd = makedist('Normal');
    pd.mu = 0;
    pd.sigma = Sigma;
    kappa = truncate(pd,-100,100);
    
    % vectors to store the averging theta values in 
    theta_bar_vec = zeros(N,1);
    theta_bar_vec_add = zeros(N,1);
    theta_bar_vec_sub = zeros(N,1);

    % vectors to store the variance of theta values in 
    theta_var_vec = zeros(N,1);
    theta_var_vec_add = zeros(N,1);
    theta_var_vec_sub = zeros(N,1);

    % draw delta_k from bernoulli +-1 dist  
    delata_k = binornd(1,0.5,1,M);
    delata_k(delata_k == 0) = -1;

    % calc c_k * delta_k 
    c_delat_k = c_m.*delata_k;
    
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

%         theta_F_log_add = theta_F(:,1) + c_delat_k(m); 
%         theta_F_log_sub = theta_F(:,1) - c_delat_k(m); 
% 
%         theta_bar_0_add =  mean(theta_F_log_add); 
%         theta_bar_0_sub =  mean(theta_F_log_sub); 


        theta_bar_0_add = mean(theta_F) + c_delat_k(m);
        theta_bar_0_sub = mean(theta_F) - c_delat_k(m);
        
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
            
%             theta_F_log_add = theta_F(:,1) + c_delat_k(m); 
%             theta_F_log_sub = theta_F(:,1) - c_delat_k(m); 
%             
%             theta_bar_vec_add(n) = mean(theta_F_log_add);
%             theta_bar_vec_sub(n) = mean(theta_F_log_sub);

            theta_bar_vec_add(n) = mean(theta_F) + c_delat_k(1);
            theta_bar_vec_sub(n) = mean(theta_F) - c_delat_k(1);
            
            % calc variance 
            if n == 1
                theta_var_vec(n) = (1/(J-1))*sum((theta_P - theta_bar_0).^2);
                theta_var_vec_add(n) = (1/(J-1))*sum((theta_P + c_delat_k(m) - theta_bar_0_add).^2);
                theta_var_vec_sub(n) = (1/(J-1))*sum((theta_P - c_delat_k(m) - theta_bar_0_sub).^2);


            else
                theta_var_vec(n) = (1/(J-1))*sum((theta_P - theta_bar_vec(n-1)).^2);
                theta_var_vec_add(n) = (1/(J-1))*sum((theta_P + c_delat_k(m) - theta_bar_vec_add(n-1)).^2);
                theta_var_vec_sub(n) = (1/(J-1))*sum((theta_P - c_delat_k(m) - theta_bar_vec_sub(n-1)).^2);

            end
        end
        
        % calc log lik
        log_lik_vec(m) = sum(log(mean(w_m,2)));
        
        % update
        
        % calc gradient 
        G_k_theta = sum( theta_var_vec.^(-1).*(theta_bar_vec(1:N) - [theta_bar_0; theta_bar_vec(1:end-1)]));
        
        % calc delta G_k 
        delta_G_k = sum( theta_var_vec_add.^(-1).*(theta_bar_vec_add(1:N) - [theta_bar_0_add; theta_bar_vec_add(1:end-1)])) - ...
                    sum( theta_var_vec_sub.^(-1).*(theta_bar_vec_sub(1:N) - [theta_bar_0_sub; theta_bar_vec_sub(1:end-1)]));
        
        % calc pre-iterated estimte of the Hessian  
        H_hat_k = 0.5*(  delta_G_k/(2*c_delat_k(m)) + transpose(  delta_G_k/(2*c_delat_k(m))  ) );
        
        % calc iterated estiamte of the Hessian
        if m == 1
            H_bar_k = (1/(m+1))*H_hat_k;
        else
            H_bar_k = (m/(m+1))*H_bar_k + (1/(m+1))*H_hat_k;
        end
        
        
        H_bar_bar_k = (H_bar_k*H_bar_k)^(1/2) + 0.1*eye(1);
        H_vec(m) = H_bar_bar_k;
        theta_vec(m+1,:) = theta_vec(m,:) - a_m(m)*((H_bar_k)^(-1))*G_k_theta;

        % update theta
        %theta_vec(m+1) = theta_vec(m) - a_m(m)*(H_bar_k^(-1))*G_k_theta;
        
        % accept update 
        if abs(theta_vec(m+1) - theta_vec(m)) > 1
            nbr_not_a = nbr_not_a +1;
            theta_vec(m+1) = theta_vec(m);
            H_bar_k = eye(1); % set hessian matrix to the identity matrix when the hessian matrix is bad 
        end
        
    end
    nbr_not_a = nbr_not_a - 1;
    theta_obs =   theta_vec(1:end-1);    
end
