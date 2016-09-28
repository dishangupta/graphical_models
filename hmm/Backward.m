% This function computes the beta table 
% and the average log-likelihood using the
% backward algorithm for a Hidden Markov Model.
function [betat, averagell] = Backward(observations, transition, emission, prior)

    T = size(observations, 2);
    N = size(emission, 1);

    betat = zeros(N,T);
    Xt = observations;

    % Backward algorithm -- Dynamic Programming
    % Initialize beta

    betat(:, T) = ones(N,1);
    count = 0;

    % Recursively compute beta at each t

    for t = T-1:-1:1

        for j = 1:N
            betat(:, t) = betat(:,t) + betat(j, t+1)*emission(j, Xt(t+1))*transition(:, j);
        end

         % Circumvent small beta values
        if (sum(betat(:, t)) < 0.001)
            betat(:, t) = betat(:, t) * 10^6;
            count = count + 1;
        end

    end

    % Compute average log-likelihood
    loglik = log(sum(prior.*emission(:, Xt(1)).*betat(:, 1))) + count*log(10^(-6));
    averagell = loglik/T;
end