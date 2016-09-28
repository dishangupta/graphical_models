% This function naively decodes the maximum likelihood  
% sequence of hidden for a Hidden Markov Model using
% the maximum likelihood state for each observation 
% independent transition probabilities.

function [ Z ] = NaiveDecoder(observations, emission)

    X = observations;
    T = size(observations, 2); % #observations
       
    Z = zeros(1,T);
    
    for t = 1:T
        [prob state] = max(emission(:, X(t)));
        Z(1, t) = state;
    end
   
        
end

