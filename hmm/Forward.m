% This function computes the alpha table 
% and the average log-likelihood using the
% forward algorithm for a Hidden Markov Model.
function [alphat, averagell] = Forward(observations, transition, emission, prior)

    T = size(observations, 2);
    N = size(emission, 1);

    alphat = zeros(N,T);
    Xt = observations;

    % Forward algorithm -- Dynamic Programming
    % Initialize alpha

    alphat(:, 1) = prior.*(emission(:, Xt(1)));
    count = 0;

    % Recursively compute alpha at each t

    for t = 2:T

        for j = 1:N
            alphat(:, t) = alphat(:, t) + (emission(:, Xt(t)).*transition(j, :)')*alphat(j, t-1);
        end

        % Circumvent small alpha values
                
        if (sum(alphat(:, t)) < 0.01)
            alphat(:, t) = alphat(:, t) * 10^6;
            count = count + 1;
        end
        
        
    end

    % Compute average log-likelihood
    loglik = log(sum(alphat(:, T))) + count*log(10^(-6));
    averagell = loglik/T;

end








