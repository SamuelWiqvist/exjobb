function [theta_m] = IF2(theta0,Y_hat, M, sigma_m,J,N)
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

    
    % print msg 
    disp('Starting IF2')
    
    % default values for M,J,N, sigma_m
    if nargin < 6
        J = 1000;
        N = length(Y_hat); 
        if nargin == 1 % the default value for M is 50 
            M = 50;        
            sigma_m = 0.1.^(0:M-1);

        end
    end
    

    % simulator for f_x0 : N(0,1), 
    % simulator for f_xn | x_n-1 ; omega : N(omgea*x_n-1,1) and
    % simulator for f yn | xn : normal(xn,2) are hard coded 
    % h_n: gaussian with mean omega_j and variance sigma_m

    
    theta_m = zeros(M,J); % pre-allocate matrix for particels 
    
    index = 1:J; % set indecies 
    
    % procidure
    for m = 1:M % loop for nbr of iterations 
        
        % initilize parameters
        if m == 1 % 
            omega_F = normrnd(theta0, sigma_m(m));
        else
            omega_F = normrnd(theta_m(m-1,:)', sigma_m(m));
        end
        
        % Initilize states
        X_F = normrnd(0,1,J,1); % set X_F
        
        omega_F_minus_1 = omega_F; % set old values
        X_F_minus_1 = X_F; % set old values
        
        for n = 1:N % loop for nbr of times steps 
           
            % Preturb parameters
            omega_P = normrnd(omega_F_minus_1, sigma_m(m)); 
            
            % Simulate prediction particles
            X_P = normrnd(omega_P.*X_F_minus_1, 1); 
            
            % Evaluate weigths
            w = normpdf(Y_hat(n), X_P , 1);
            w = w./sum(w); % normalize 
            
            % Systematic resampling of indecies
            k = datasample(index,J,'Replace',true, 'Weights', w); 
            
            % Resaple 
            omega_F = omega_P(k);
            X_F = X_P(k);
            
            % Set old values to current values
            omega_F_minus_1 = omega_F; 
            X_F_minus_1 = X_F;
        end
        % Store new parameter swarm
        theta_m(m,:) = omega_F; 
    end
    
    % Add the inital parameter swarm
    theta_m = [theta0'; theta_m]; 
end