function [theta_obs, log_lik_vec, nbr_not_a, eig_vec, norm_vec] = IF1_gompertz_v2_Hessian(theta0_log,Y_hat, M, N,J_seq,tau_seq,sigma_seq,Sigma,par_K, print_m, a_m,c_m,burn_in)
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
    disp('Starting IF1_gompertz_v2_A')
        
    eig_vec = zeros(3,M);
    norm_vec = zeros(M,1);
    % set start value for the nbr of steps that are not accepted
    nbr_not_a = 1;
    
    log_lik_vec = zeros(M,1);
    
    theta_vec = zeros(M+1,3);
    theta_vec(1,:) = theta0_log;
    
    % set distribution for parameter perturbation 
    pd_r = makedist('Normal');
    limit = 1000;
    pd_r.mu = 0;
    pd_r.sigma = Sigma(1,1);
    kappa_r = truncate(pd_r,-limit,limit);
    
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
    
    % vectors to store the averging theta values in 
    theta_bar_vec_log = zeros(N,3);
    theta_bar_vec_add_log = zeros(N,3);
    theta_bar_vec_sub_log = zeros(N,3);

    % vectors to store the variance of theta values in 
    theta_var_vec_log = zeros(N,3);
    theta_var_vec_add_log = zeros(N,3);
    theta_var_vec_sub_log = zeros(N,3);

    % draw delta_k from bernoulli +-1 dist  
    delata_k_r = binornd(1,0.5,1,M);
    delata_k_sigma = binornd(1,0.5,1,M);
    delata_k_theta = binornd(1,0.5,1,M);
    delata_k = [delata_k_r', delata_k_sigma', delata_k_theta'];   
    delata_k(delata_k == 0) = -1;

    % calc c_k * delta_k 
    c_delat_k = [c_m'.*delata_k(:,1) c_m'.*delata_k(:,2) c_m'.*delata_k(:,3)];
    
    

    
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
        Z_r = random(kappa_r,J,1);
        Z_sigma = random(kappa_sigma,J,1);
        Z_tau = random(kappa_tau,J,1);
        
        % Initilize parameters
        tau_Z = tau_seq(m).*[Z_r, Z_sigma, Z_tau];
        theta_F_log = [theta_vec(m,1) + tau_Z(:,1), theta_vec(m,2) + tau_Z(:,2), theta_vec(m,3) + tau_Z(:,3)];
        
        % Initilize states 
        X_F = rand(J,1);
        
        % calc theta_bar_0
        theta_bar_0_log = mean(theta_F_log);
        
%         theta_F_log_add = [theta_F_log(:,1) + c_delat_k(m,1), theta_F_log(:,2) + c_delat_k(m,2), theta_F_log(:,3) + c_delat_k(m,3)]; 
%         theta_F_log_sub = [theta_F_log(:,1) - c_delat_k(m,1), theta_F_log(:,2) - c_delat_k(m,2), theta_F_log(:,3) - c_delat_k(m,3)]; 
%         
%         theta_bar_0_add_log = mean(theta_F_log_add);
%         theta_bar_0_sub_log = mean(theta_F_log_sub);
        
        theta_bar_0_add_log = mean(theta_F_log) + c_delat_k(m,:);
        theta_bar_0_sub_log = mean(theta_F_log) - c_delat_k(m,:);
        
        for n = 1:N
            % Perturb parameters
            Z_r = random(kappa_r,J,1);
            Z_sigma = random(kappa_sigma,J,1);
            Z_tau = random(kappa_tau,J,1);
            
            theta_P_log = theta_F_log + sigma_seq(m).*[Z_r, Z_sigma, Z_tau];
            
            % Simulate prdiction particles
            eps_temp = exp(normrnd(0,real(exp(theta_P_log(:,2)))));
            X_P = par_K.^(1-exp(-real(exp(theta_P_log(:,1))))).*X_F.*(exp(-real(exp(theta_P_log(:,1))))).*eps_temp;%normrnd(theta_P.*X_F,1);
            
            % Evaluate weigths
            w = normpdf(log(Y_hat(n)), log(X_P) , real(exp(theta_P_log(:,3)))); 
            prob_resampling = w./sum(w); % normalize calc resampling probability 
            
            % store wegiths
            
            w_m(n,:) = w; 
            % Systematic resampling of indecies 
            k = datasample(index,J,'Replace',true, 'Weights', prob_resampling);
            
            % Resaple 
            X_F = X_P(k);
            theta_F_log = theta_P_log(k,:);
            
            % calc mean
            theta_bar_vec_log(n,:) = mean(theta_F_log);
%             
%             theta_F_log_add = [theta_F_log(:,1) + c_delat_k(m,1), theta_F_log(:,2) + c_delat_k(m,2), theta_F_log(:,3) + c_delat_k(m,3)]; 
%             theta_F_log_sub = [theta_F_log(:,1) - c_delat_k(m,1), theta_F_log(:,2) - c_delat_k(m,2), theta_F_log(:,3) - c_delat_k(m,3)]; 
%             
%             theta_bar_vec_add_log(n,:) = mean(theta_F_log_add);
%             theta_bar_vec_sub_log(n,:) = mean(theta_F_log_sub);
            
            theta_bar_vec_add_log(n,:) = mean(theta_F_log) + c_delat_k(m,:);
            theta_bar_vec_sub_log(n,:) = mean(theta_F_log) - c_delat_k(m,:);
            
            
            % calc variance 
            if n == 1
                theta_var_vec_log(n,:) = [(1/(J-1))*sum((theta_P_log(:,1) - theta_bar_0_log(1)).^2), ...
                                          (1/(J-1))*sum((theta_P_log(:,2) - theta_bar_0_log(2)).^2) ,...
                                          (1/(J-1))*sum((theta_P_log(:,3) - theta_bar_0_log(3)).^2)];                                  
                
                theta_var_vec_add_log(n,:) = [(1/(J-1))*sum((theta_P_log(:,1) + c_delat_k(m,1) - theta_bar_0_add_log(1)).^2), ....
                                              (1/(J-1))*sum((theta_P_log(:,2) + c_delat_k(m,2) - theta_bar_0_add_log(2)).^2) ,...
                                              (1/(J-1))*sum((theta_P_log(:,3) + c_delat_k(m,3) - theta_bar_0_add_log(3)).^2)];
                
                
                theta_var_vec_sub_log(n,:) = [(1/(J-1))*sum((theta_P_log(:,1) - c_delat_k(m,1) - theta_bar_0_sub_log(1)).^2), ....
                                              (1/(J-1))*sum((theta_P_log(:,2) - c_delat_k(m,2) - theta_bar_0_sub_log(2)).^2) ,...
                                              (1/(J-1))*sum((theta_P_log(:,3) - c_delat_k(m,3) - theta_bar_0_sub_log(3)).^2)];
                                        

            else     
                theta_var_vec_log(n,:) = [(1/(J-1))*sum((theta_P_log(:,1) - theta_bar_vec_log(n-1,1)).^2), ...
                                          (1/(J-1))*sum((theta_P_log(:,2) - theta_bar_vec_log(n-1,2)).^2), ...
                                          (1/(J-1))*sum((theta_P_log(:,3) - theta_bar_vec_log(n-1,3)).^2)];
                                  
                
                theta_var_vec_add_log(n,:) = [(1/(J-1))*sum((theta_P_log(:,1) + c_delat_k(m,1) - theta_bar_vec_add_log(n-1,1)).^2), ....
                                              (1/(J-1))*sum((theta_P_log(:,2) + c_delat_k(m,2) - theta_bar_vec_add_log(n-1,2)).^2) ,...
                                              (1/(J-1))*sum((theta_P_log(:,3) + c_delat_k(m,3) - theta_bar_vec_add_log(n-1,3)).^2)];
                
                
                theta_var_vec_sub_log(n,:) = [(1/(J-1))*sum((theta_P_log(:,1) - c_delat_k(m,1) - theta_bar_vec_sub_log(n-1,1)).^2), ....
                                              (1/(J-1))*sum((theta_P_log(:,2) - c_delat_k(m,2) - theta_bar_vec_sub_log(n-1,2)).^2) ,...
                                              (1/(J-1))*sum((theta_P_log(:,3) - c_delat_k(m,3) - theta_bar_vec_sub_log(n-1,3)).^2)];
                                        
            end
        end
        
        % calc log lik
        log_lik_vec(m) = sum(log(mean(w_m,2)));
        
        % update

        % calc gradient 
        theta_var_vec_log = theta_var_vec_log.^(-1);
        
        G_k_theta = [sum( theta_var_vec_log(:,1).*(theta_bar_vec_log(1:N,1) - [theta_bar_0_log(1); theta_bar_vec_log(1:end-1,1)])), ...
                     sum( theta_var_vec_log(:,2).*(theta_bar_vec_log(1:N,2) - [theta_bar_0_log(2); theta_bar_vec_log(1:end-1,2)])), ...
                     sum( theta_var_vec_log(:,3).*(theta_bar_vec_log(1:N,3) - [theta_bar_0_log(3); theta_bar_vec_log(1:end-1,3)]))]';
          
        % calc delta G_k 
        theta_var_vec_add_log_r = theta_var_vec_add_log(:,1).^(-1);
        theta_var_vec_add_log_sigma = theta_var_vec_add_log(:,2).^(-1);
        theta_var_vec_add_log_tau = theta_var_vec_add_log(:,3).^(-1);

        theta_var_vec_sub_log_r = theta_var_vec_sub_log(:,1).^(-1);
        theta_var_vec_sub_log_sigma = theta_var_vec_sub_log(:,2).^(-1);
        theta_var_vec_sub_log_tau = theta_var_vec_sub_log(:,3).^(-1);

        delta_G_k = [sum( theta_var_vec_add_log_r.*(theta_bar_vec_log(1:N,1) - [ theta_bar_0_add_log(1); theta_bar_vec_log(1:end-1,1)])) - ...
                    sum( theta_var_vec_sub_log_r.*(theta_bar_vec_log(1:N,1) - [ theta_bar_0_sub_log(1); theta_bar_vec_log(1:end-1,1)])), ...
                    sum( theta_var_vec_add_log_sigma.*(theta_bar_vec_log(1:N,2) - [theta_bar_0_add_log(2); theta_bar_vec_log(1:end-1,2)])) - ...
                    sum( theta_var_vec_sub_log_sigma.*(theta_bar_vec_log(1:N,2) - [theta_bar_0_sub_log(2); theta_bar_vec_log(1:end-1,2)])), ...
                    sum( theta_var_vec_add_log_tau.*(theta_bar_vec_log(1:N,3) - [theta_bar_0_add_log(3); theta_bar_vec_log(1:end-1,3)])) - ...
                    sum( theta_var_vec_sub_log_tau.*(theta_bar_vec_log(1:N,3) - [theta_bar_0_sub_log(3); theta_bar_vec_log(1:end-1,3)]))]';
                
        

        
        % calc pre-iterated estimte of the Hessian  
        %H_hat_k = 0.5*(  transpose(delta_G_k')/(2*c_delat_k(m,:)') + transpose(  transpose(delta_G_k)/(2*c_delat_k(m,:)')  ) );
        
        H_hat_k = 0.5*(  delta_G_k*((2*c_delat_k(m,:)).^(-1)) + transpose(  delta_G_k*((2*c_delat_k(m,:)).^(-1)) ) );
        
        % calc iterated estiamte of the Hessian
        if m == 1
            H_bar_k = (1/(m+1))*H_hat_k;
        else
            H_bar_k = (m/(m+1))*H_bar_k + (1/(m+1))*H_hat_k;
        end
        
        d_k_vec = 10000*(1:1:M);
        
        d_k = d_k_vec(m);
        
        %H_bar_bar_k = real(sqrtm(H_bar_k*H_bar_k)) + d_k*eye(3);
        H_bar_bar_k = H_bar_k + d_k*eye(3);
        
        eig_vec(:,m) = eig(H_bar_bar_k);
        norm_vec(m) = norm(H_bar_bar_k, 'fro');
        
        % update theta
        if m >= burn_in
            theta_vec(m+1,:) = theta_vec(m,:) - a_m(m)*(H_bar_bar_k\G_k_theta)';
        else
             theta_vec(m+1,:) = theta_vec(m,:) + a_m(m)*(G_k_theta)';
        end
        
        % accept update
        if abs(theta_vec(m+1,1) - theta_vec(m,1)) > 1 ||  abs(theta_vec(m+1,2) - theta_vec(m,2)) > 1 ||  abs(theta_vec(m+1,2) - theta_vec(m,2)) > 1
            theta_vec(m+1,:) = theta_vec(m,:);
            %H_bar_k = eye(3);
            nbr_not_a = nbr_not_a + 1;
        end
        
        
        
    end
    nbr_not_a = nbr_not_a - 1;
    theta_obs =   theta_vec(1:end-1,:);    
end
