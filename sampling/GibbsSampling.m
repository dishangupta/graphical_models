% Forward Sampling -- Initialize parameters
mu = [-5 5];
var = [1 1];
tau = 10;
p = 1/2;
N = 100;

% Sample x
x = GaussianMixture(mu, var, p, N);

% Gibbs Sampling -- Initialize parameters
z = size(x);
mu = [-1 1];
burn = 500;
numSamples = 2000;
MU = zeros(numSamples-burn, size(mu, 2));

% Generate 2000 samples

for s = 1:numSamples
    % Sample Z's from mu's
    for i = 1:N
        pZ = normpdf(x(i), mu, var);
        pZ = pZ./(sum(pZ));
        z(i) = Bernoulli(pZ(1));
    end

    % Sample mu's from Z's
    for j = 1:size(mu, 2)
        xj = x(find(z == j));
        nj = numel(find(z == j));
        postMu = sum(xj)/(nj + 1/tau);
        postVar = 1/(nj + 1/tau);
        mu(j) = mvnrnd(postMu, postVar);
    end
    
    % Burn 500 samples
    if (s > burn)
        MU(s-burn, :) = mu;
    end
    
end

% Plot kernel density based on a normal smooth kernel
[f1, xi1] = ksdensity(MU(:,1));
[f2, xi2] = ksdensity(MU(:,2));
plot(xi1, f1, xi2, f2); 



 