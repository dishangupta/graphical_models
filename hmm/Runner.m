% Main Script to train, decode and evaluate
% and HMM

N = 2; % #states
O = 27; % #output symbols

fileTrain = 'hmm-train-cleaned.txt';
fileDecode = 'hmm-decode-cleaned.txt';
fileTest = 'hmm-test-cleaned.txt';

% Train HMM
disp(['============================================HMM Training============================================' char(10) char(10)]);
X = TextToNumeric(fileTrain); % observations
numRuns = 20;
convergence = 0.00001;
[transition, emission, prior] = HMMTrain(X, N, O, numRuns, convergence);

% Decode HMM
disp(['============================================HMM Decoding============================================' char(10)]);
X = TextToNumeric(fileDecode); % observations
Z = Viterbi(X, transition, emission, prior);
disp(['The maximum likelihood hidden state sequence using the Viterbi algorithm is:' char(10) num2str(Z) char(10) char(10)]);

% Test HMM
disp(['============================================HMM Testing============================================' char(10) char(10)]);
X = TextToNumeric(fileTest); % observations
[alpha averagell] = Forward(X, transition, emission, prior);
display(['The average log-likelihood for the test set is: ' num2str(averagell) char(10)]);








