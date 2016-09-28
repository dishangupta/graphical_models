% This function updates the HMM parameters prior,
% transition and emission using the maximization step of the
% Baum-Welch algorithm for training a Hidden Markov 
% Model.

function [prior, transition, emission] = Maximization(eta, gamma, observations, numSymbols )
    
    X = observations;
    N = size(eta, 1);
    T = size(X, 2);
    O = numSymbols;
    
    % Update prior
    prior = gamma(:, 1);
    prior = prior./repmat(sum(prior, 1), N, 1);
    
    % Update transition matrix
    transition = sum(eta, 3)./repmat(sum(gamma(:, 1:T-1), 2), 1, N);
    transition = transition./repmat(sum(transition, 2), 1, N); % Normalize 
    
    % Update emission matrix
    emission = zeros(N, O);
    for o = 1:O
        ind = find(X == o); 
        % Smooth for zero count symbols in training
        emission(:, o) = (sum(gamma(:, ind), 2) + 1)./(sum(gamma, 2) + 1); 
    end
    emission = emission./repmat(sum(emission, 2), 1, O); % Normalize

end

