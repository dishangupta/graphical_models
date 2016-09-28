% Function to run the neigherborhood estimation algorithm:
% f(theta_a) = 0.5*||Xa - X/a*theta_a||_2^2 + lambda*||theta_a||_1

% Inputs:
%   X: n-by-p input matrix of n p-dimensional data points
%   lambda: weight for l1 norm penalty

% Outputs:
%   Theta: estimate (sparse) of the inverse covariance matrix

function [ Theta ] = Neigbhorhood_Estimate( X, lambda )

    %  Intialize values
    n = size(X,1);
    p = size(X,2); 
    Theta = zeros(p, p);

    % Normalize data to zero-mean distribution
    mu = sum(X,1)/n;
    X = X - repmat(mu, n, 1);

    for a = 1:p
        Xa = X(:, a);
        I = ~ismember(1:p, a);
        Xna = X(:, I);
        theta_a = lasso(Xna, Xa, lambda);
        theta_a = [theta_a(1:a-1); 0; theta_a(a:end)];
        Theta(a, :) = theta_a; 
    end

    % Symmetrize Theta using OR rule
    for i = 1:p
        for j = 1:p
            if(Theta(i,j) == 0)
                Theta(j,i) = 0;
            end
        end
    end
    
   % Plot binary image 
   Plot(Theta);
end