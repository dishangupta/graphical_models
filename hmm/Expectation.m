% This function computes the temporary variables
% eta and gamma using the expectation step of the
% Baum-Welch algorithm for training a Hidden Markov 
% Model.

function [eta, gamma] = Expectation(observations, alpha, beta, transition, emission)
    
    T = size(observations, 2); % #observations
    N = size(transition, 1); % states  
    X = observations;

    % Initialize temporary variables
    eta = zeros(N, N, T-1); 
    gamma = zeros(N, T);
    
    for t = 1:T-1
        den = 0;
        for i = 1:N
            for j = 1:N
                eta (i, j, t) = alpha(i, t)*transition(i,j)*emission(j, X(t+1))*beta(j, t+1);
            end
            den =  den + alpha(i, t)*beta(i, t);
            gamma(i, t) = alpha(i, t)*beta(i, t);
        end

        eta(:, :, t) = eta(:, :, t)/den;
        gamma(:, t) = gamma(:, t)/den;
    end
    
     % Update emission matrix
    gamma(1, T) = alpha(1, T)*beta(1, T)/(alpha(1, T)*beta(1, T) + alpha(2, T)*beta(2, T));
    gamma(2, T) = alpha(2, T)*beta(2, T)/(alpha(1, T)*beta(1, T) + alpha(2, T)*beta(2, T));
end

