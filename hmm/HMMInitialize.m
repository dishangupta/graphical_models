% This function randomly initializes the prior, transition
% and emission parameters of a Hidden Markov Model.

function [transition, emission, prior] = HMMInitialize(numStates, numSymbols)

    N = numStates;
    O = numSymbols;

    % Randomly initialize HMM parameters 
    transition = rand(N, N);
    emission = rand(N, O);
    prior = rand(N, 1);

    % Normalize Initializations
    transition = transition./repmat(sum(transition, 2), 1, N);
    emission = emission./repmat(sum(emission, 2), 1, O);
    prior = prior./repmat(sum(prior, 1), N, 1);
    

end

