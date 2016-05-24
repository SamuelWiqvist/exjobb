function [theta_hat_ML] = ML_est(Y_hat,theta0)
    % Maximumlikelihood estimation for the AR(1) model. 
    % Returns the ML estimations of \theta_1 
    disp('Starting ML')
    %theta0 = -1 + (1+1)*rand(1,1);
    % use the matlab built in function to estimate the parameters using ML'
    % specify model
    A = NaN;
    B = 1;
    C = 1;
    D = 1;
    Mdl = ssm(A,B,C,D);
    
    [EstMdl,EstParamCov,logL,info] = estimate(Mdl,Y_hat, theta0);
    theta_hat_ML = EstMdl.A;
end