% Function to run the Glasso algorithm:
% f(beta) = 0.5*||W11^(1/2)*beta - W11^(-1/2)*s12||_2^2 + lambda*||beta||_1

% Inputs:
%   X: n-by-p input matrix of n p-dimensional data points
%   lambda: weight for l1 norm penalty

% Outputs:
%   Theta: estimate (sparse) of the inverse covariance matrix


function [ Theta ] = Glasso( X, lambda )

    %  Dimensions
    n = size(X,1);
    p = size(X,2); 

    % Normalize data to zero-mean distribution
    mu = sum(X,1)/n;
    X = X - repmat(mu, n, 1);

    % Initialize values
    S = (1/n)*((X')*(X));
    W = S + lambda*eye(p);
    Wnew = zeros(size(W));
    Theta = zeros(p, p);
    Theta_new = zeros(p,p);

    % set options for optimization
    miniter = 10;
    maxiter = 5000;
    tol = 1e-5;

    % Implement Glasso
    for iter = 1:maxiter
        for a = 1:p
            I = ~ismember(1:p, a);
            W11 = W(I, I);
            w22 = W(a,a);
            w12 = W(a,:);
            w12 = (w12(I))';
            S11 = S(I,I);
            s22 = S(a,a);
            s12 = S(a,:);
            s12 = (s12(I))';

            % Compute W11^(1/2)
            [V,D] = eig(W11);
            W11_sqrt = V * sqrt(D) * V';
            W11_sqrt_inv = inv(W11_sqrt);

            % Execute lasso
            b = (W11_sqrt_inv)*s12;
            beta = lasso(W11_sqrt, b, lambda);

            % Computer row/column of Wnew and Beta
            w12 = W11*beta;
            col_a = [w12(1:a-1); w22; w12(a:end)];
            Wnew(a,:) = col_a';
            Wnew(:,a) = col_a;
           
            theta22 = 1/(w22 - (w12')*beta);
            theta12 = -beta*theta22;
            
            Theta_new(a,:) = [theta12(1:a-1); theta22; theta12(a:end)]'; 
            Theta_new(:,a) = [theta12(1:a-1); theta22; theta12(a:end)]; 
        end

        diff = Theta_new - Theta;
        if (iter >= miniter && abs(sum(diff(:))) < tol)
            break;
        end
        
        W = Wnew;
        Theta = Theta_new;
    end
    
    % Plot binary image
    Plot(Theta);
end

