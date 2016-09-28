% This function decodes the maximum likelihood  
% sequence of hidden for a Hidden Markov Model using
% the Viterbi algorithm.

function [Z] = Viterbi(observations, transition, emission, prior)
    
    X = observations;
    T = size(observations, 2); % #observations
    N = size(transition, 1); % #states
    
    maxPath = zeros(N,T);
    maxValues = zeros(N,T);
    Z = zeros(1,T);
    
    % Viterbi algorithm -- Dynamic Programming
    % Initialize mu's
    mu = prior.*emission(:, X(1)); 

    for t = 2:T

        for j = 1:N
            objective = emission(j, X(t))*(transition(:, j).*(mu));
            [value, ind] = max(objective);
            maxValues(j, t-1) = value;
            maxPath(j, t-1) = ind;
        end

        mu = maxValues(:, t-1);

    end

    [value, ind] = max(mu);
    Z(T) = ind;

    % Backtrace
    for t = T-1:-1:1
        Z(t) = maxPath(Z(t+1), t);
    end
  
end


