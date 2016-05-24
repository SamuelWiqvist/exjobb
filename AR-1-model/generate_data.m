function [X,Y] = generate_data(N)
    %Linear POMP model for IF1 and IF2, same exaple as in (Inoides 2011)
    x0 = 0;
    X = -1 + (1+1)*rand(N,1);
    Y = zeros(N,1);
    eps = normrnd(0,1,N,1);
    eta = normrnd(0,1,N,1);
    omega = 0.8; % true parameter

    % generate data
    X(1) = omega*x0 + eps(1);
    Y(1) = X(1) + eta(1);
    for i = 2:N    
        X(i) = omega*X(i-1) + eps(i);
        Y(i) = X(i) + eta(i);
    end
    %Plot data 
    %figure
    %plot(X, 'b')
    %hold on
    %plot(Y, 'r')
end
