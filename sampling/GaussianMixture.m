%Sample from mixture of Gaussians

function [x] = GaussianMixture(mu, var, p, N)
   
    x = size(N,1);
    
    for i = 1:N 
        % Sample Z ~ Bernoulli(p)
        z = Bernoulli(p);

        % Sample values (Xi = xi) based on samples of Z
        if (z == 1)
            x(i) = mvnrnd(mu(1), var(1));
        else 
            x(i) = mvnrnd(mu(2), var(2));
        end
    end
    
end

