% Load Data
load('data.mat');

% Initialize
M = size(data, 1);
K = size(beta_matrix, 2);
Theta = zeros(M, K);
alpha = 0.01;
totalIter = 0;

tic
for m = 1:M
    
    [phi_ret gamma_ret, iter_ret] = IndividualVariationalInference(data, beta_matrix, alpha, m);
    
    Theta(m, :) = gamma_ret';
    totalIter = totalIter + iter_ret;

end
toc

imagesc(Theta);
colorbar;

%imagesc(gamma_ret);
%colorbar;



