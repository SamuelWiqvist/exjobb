function [omega_hat_vec, hessian_approx_vec] = ExtendedKalmanFilter(omega0, Y_hat)
    % Estimates the parameters using the EKF algorithm 
    
    disp('Starting EKF')

    x0 = -1 + (1+1)*rand(1,1); 
    N = length(Y_hat);

    syms x omega

    f(x, omega) = [x*omega, omega]';
    F = jacobian(f); 

    h(x,omega) = [x]';                   
    H = jacobian(h);

    f = matlabFunction(f);
    F = matlabFunction(F);
    h = matlabFunction(h);
    H = matlabFunction(H);


    yt = Y_hat; % Observations
    xt = zeros(N,2); % Our HMM
    Pt = zeros(N,2); 
    
    hessian_approx_vec = zeros(N,1);
    
    % Initial guess
    xtt = [x0, omega0]'; %  %x1,x2, k,s,a,b,c
                                               % A good set of parameters
    xt(1,:) = xtt;

    R = 0.001^2*eye(1); % R = Rw, Q = Re and is time dependent
    Pxx = 0.005*eye(2); % Rxx

    for i = 2:N
        % Time dependent variables 
        q = 0.001^2;
        Q =  diag([1, q]); % Re

        % Predict
        xtt_1 = f(xtt(1), xtt(2));
        Ft = F(xtt(1), xtt(2));
        Pxx_1 = Ft*Pxx*Ft' + Q; % Rxx_1

        % Update
        y_tilde = yt(i,:)' - h(xtt_1(1), xtt_1(2));
        Ht = H(xtt_1(1), xtt(2));
        St = Ht * Pxx_1 * Ht' + R; % Ryy
        Kt = (Pxx_1 * Ht') / St; % Rxx
        Pxx = Pxx_1 - Kt*Ht*Pxx_1;
        xtt = xtt_1 + Kt*y_tilde;
        
        %save
        hessian_approx_vec(i) = St;
        xt(i,:) = xtt;
        Pt(i,:) = diag(Pxx);
    end
    omega_hat_vec = xt(:,2);
end
