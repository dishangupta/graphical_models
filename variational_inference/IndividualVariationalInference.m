function [ phi_ret, gamma_ret, iter_ret ] = IndividualVariationalInference( data, beta_matrix, alph, doc )

% Run Variational Inference for 1 individual (document)

    % Initialize gamma and phi
    N = size(data, 2); 
    K = size(beta_matrix, 2);
    ind = find(data(doc,:) ~= 0);
    n1 = numel(ind);
    alpha = alph*ones(K, 1);
    phi = (1/K)*ones(n1, K);
    gamma = alpha + N/K ;

    % Run iterations till convergence
    eps = 0.001;
    conv = 1000;
    iter = 0;

    phi_temp = zeros(size(phi));
    phi_new = zeros(size(phi));

    while (1)
        iter = iter + 1;
        
        % Update phi
        for n = 1:n1
            for k = 1:K
                phi_new(n, k) = beta_matrix(ind(n), k)*exp(psi(gamma(k)));  
            end
        end

        % Normalize phi
        for n = 1:n1
            phi_new(n, :) = phi_new(n, :)./sum(phi_new(n, :));
        end

        for k = 1:K
            phi_temp(:, k) = phi_new(:,k).*data(doc, ind)';
        end

        % Update gamma
        gamma_new = alpha + (sum(phi_temp, 1))';


        % Check convergence
        a = abs(gamma_new - gamma);
        b = abs(phi_new - phi);

        if (all(a < eps) & all(b < eps))
            break;
        end

        phi = phi_new;
        gamma = gamma_new;

    end
    
    % Return values
    phi_ret = phi_new;
    gamma_ret = gamma_new;
    iter_ret = iter;
end

