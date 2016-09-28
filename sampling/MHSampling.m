% Forward Sampling -- Initialize parameters
mu = [-5 5];
var = [1 1];
p = 1/2;
N = 100;

% Sample x
x = GaussianMixture(mu, var, p, N);

% Metropolis-Hastings -- Initialize parameters
tau_prior = 10;
burn = 500;
numSamples = 2000;
MU = zeros(numSamples-burn, size(mu, 2));

mu = [-1 1];
Tau = zeros(size(mu,2), size(mu,2)); 
for i = 1:size(Tau, 1)
    for j = 1:size(Tau, 2)
        if (i == j)
            Tau(i,j) = 10;
        else
            Tau(i,j) = 4;
        end    
    end
end

% Assume proposal distribution Q = N(mu'|mu, Tau) and 
% generate sample from it
s = 0;
while (s ~= 2000)
    
    mu_next = mvnrnd(mu, Tau);

    % Compute acceptance probability A(mu -> mu')
    A = 1;
       
    num1 = normpdf(x, mu_next(1), var(1));
    num2 = normpdf(x, mu_next(2), var(2));
    den1 = normpdf(x, mu(1), var(1));
    den2 = normpdf(x, mu(2), var(2));
    
    A = (p*num1 + (1-p)*num2)./ (p*den1 + (1-p)*den2);
    A = prod(A);
   
    % Incorporate prior
    ratio = (normpdf(mu_next(1), 0, tau_prior)*normpdf(mu_next(2), 0, tau_prior))/(normpdf(mu(1), 0, tau_prior)*normpdf(mu(2), 0, tau_prior));
    A = A * ratio;

    if (A > 1)
        A = 1;
    end
    
    % For acceptance
    if (~isnan(A) && Bernoulli(A) == 1)
        % Burn 500 samples
        s = s+1;
        
        if (s > burn)
            MU(s-burn, :) = mu_next;
        end
        
        mu = mu_next;
    end
    
end

% Plot kernel density based on a normal smooth kernel
[f1, xi1] = ksdensity(MU(:,1));
[f2, xi2] = ksdensity(MU(:,2));
plot(xi1, f1, xi2, f2); 
