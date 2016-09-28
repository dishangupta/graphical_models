% This function implements the Baum-Welch Expectation-
% Maximization algorithm to train the prior, transition
% and emission parameters of a Hidden Markov Model.

function [transition, emission, prior] = HMMTrain(observations, numStates, numSymbols, numRuns, convergence)

    X = observations;
    N = numStates;
    O = numSymbols;
    convergedAll = -Inf;
    
    % Run Baum-Welch for numRuns iterations to find best parameters 
    for run = 1:numRuns    
        
        disp(['---------------------------Run' num2str(run) '---------------------------' char(10)]);
        disp(['Printing average log-likelihood at each iteration:' char(10)]);
        
        [transition, emission, prior] = HMMInitialize(N, O);

        [alpha averagell] = Forward(X, transition, emission, prior);
        [beta] = Backward(X, transition, emission, prior);

        disp([char(9) num2str(averagell)]);
                
        % Iterate till convergence
        while (1)

            % Compute temporary variables    
            [eta gamma] = Expectation(X, alpha, beta, transition, emission);

            % Update HMM parameters
            [prior transition emission] = Maximization(eta, gamma, X, O);

            [alpha averagellNew] = Forward(X, transition, emission, prior);
            [beta] = Backward(X, transition, emission, prior);
    
            % Check convergence of EM
            if (abs(averagellNew - averagell) < convergence)
                break;
            end

            averagell = averagellNew;
            disp([char(9) num2str(averagell)]);
        end
        
        convergedAllNew = averagell;
        
        % Check for higher local maxima
        if(convergedAllNew > convergedAll)
           s = struct('transition', transition, 'emission', emission, 'prior', prior);
           convergedAll = convergedAllNew;
        end
        
        disp([char(10) char(10)]);
           
    end
   
   transition = s.transition;
   emission = s.emission;
   prior = s.prior;
   disp(['The highest average log-likelihood at convergence is: ' num2str(convergedAll)]);
end

