function [theta_m] = IF2_gompertz_est_var(theta0_log,Y_hat, M, sigma_m,J,N,par_K)
    % Description
    % Estimates the parameters theta using IF2 agirithm described in 
    % Ionides (2015)
    %
    % Inputs: 
    % theta0 = inital particle swarm
    % y_hat = obseravtions 
    % M = nbr of iterations
    % sigma_m: preturbation sequance (decreses exponentially)
    % J: nbr of particles
    % N: nbr of time steps 

    Sigma = eye(2);
    % print msg 
    disp('Starting IF2')
    
    r = 0.1;

    % simulator for f_x0 : N(0,1), 
    % simulator for f_xn | x_n-1 ; omega : N(omgea*x_n-1,1) and
    % simulator for f yn | xn : normal(xn,2) are hard coded 
    % h_n: gaussian with mean omega_j and variance sigma_m

    
    theta_m = zeros(M,J*2); % pre-allocate matrix for particels 
    
    index = 1:J; % set indecies 
    
    % procidure
    for m = 1:M % loop for nbr of iterations 
        % initilize parameters
        if m == 1 % 
            theta_F = mvnrnd(theta0_log, sigma_m(m)*Sigma);
        else
            theta_F = mvnrnd([theta_m(m-1,(1:J)) ; theta_m(m-1,J+1:2*J) ]', sigma_m(m)*Sigma);
        end
        
        % Initilize states
        X_F = rand(J,1); % set X_F
                
        for n = 1:N % loop for nbr of times steps 
           
            % Preturb parameters
            theta_P = mvnrnd(theta_F, sigma_m(m)*Sigma); 
            
            % Simulate prediction particles
            %X_P = normrnd(exp(theta_P_log).*X_F_minus_1, 1); 
            eps_temp = exp(normrnd(0,real(exp(theta_P(:,1)))));
            X_P = par_K.^(1-exp(-r)).*X_F.*(exp(-r)).*eps_temp;%normrnd(theta_P.*X_F,1);

            % Evaluate weigths
            w = normpdf(log(Y_hat(n)), log(X_P) , real(exp(theta_P(:,2))));
            prob_resampling = w./sum(w); % normalize 
            
            % Systematic resampling of indecies
            k = datasample(index,J,'Replace',true, 'Weights', prob_resampling); 
            
            % Resaple 
            theta_F = theta_P(k,:);
            X_F = X_P(k,:);
        end
        % Store new parameter swarm
        theta_m(m,:) = [theta_F(:,1)' theta_F(:,2)'] ;
    end
    
    % Add the inital parameter swarm
    % theta_m = [theta0_log'; theta_m]; 
end