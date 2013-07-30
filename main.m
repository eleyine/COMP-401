%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% General Parameters
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

WORKINGDIR = '/Users/Eleyine/Dropbox/Winter2013/COMP401/code/project';
DATAFILENAME = '/Users/Eleyine/Dropbox/Winter2013/COMP401/code/data/SBT_ETTCPAP_data.mat';
SUBSAMPLE_RATIO = 0.001;
DISCRETISE = 0;
FEATURES = [
    'fmx'; % Instanteneous respiratory frequency (Hz) sampled at 50 Hz
    'mab'; % Movement artifact statistic for Abdomen sampled at 50 Hz
    'mrc'; % Movement artifact statistic for Ribcage sampled at 50 Hz
%     'rms'; % Root mean square statistic sampled at 50 Hz
%     'prc'; % Pause statistic for Ribcage sampled at 50 Hz
%     'pab'; % Pause statistic for Abdomen sampled at 50 Hz
%     'phi'; % Asynchronic statistic sampled at 50 Hz
%     'fxp'; % Power on fmx sampled at 50 Hz
%     'fxv'; % Variance of fmx sampled at 50 Hz
%     'hmx'; % Heart rate estimate sampled at 50 Hz
];

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Initial workspace setup
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
fprintf('Loading Data...\n');

addpath(genpath(WORKINGDIR));

% avoid loading data if it is already in workspace
if ~any(strcmp(who,'Features')) || ~exist('Features')
    load(DATAFILENAME, 'Features', 'Index*', 'Fs');
end

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Methods outlined in report, uncomment desired method
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% METHOD 1 %
%----------%
PARAMS = struct('timeFrame', [-3 -2], ...
                'samplingTime', 2, ...
                'workspace', 'method1.mat', ...
                'maxIter', 10, ...
                'M', 2, ...
                'Q', 3, ...
                'method', 'forest'); % democracy %

% METHOD 2
% PARAMS = struct('timeFrame', [-23 -1]);

% COMMON %
%--------%
PARAMS.noNaN = 1; % remove trailing and leading NaN
PARAMS.ignoreSamples = 5:55; % baby samples to be ignored
PARAMS.indexPos = IndexSucced(~ismember(IndexSucced, PARAMS.ignoreSamples));
PARAMS.indexNeg = IndexFailed(~ismember(IndexFailed, PARAMS.ignoreSamples));
PARAMS.sf = Fs;
features = getSubsetStruct(Features, FEATURES);

fprintf('Subsampling Data...\n');
% Method call from 'subsample.m'
data = subsample(features, PARAMS);

fprintf('Cross-validation...\n');

% Partition indices into 5
allIndex = vertcat(PARAMS.indexPos, PARAMS.indexNeg);
cv = cvpartition(allIndex, 'k', 5)
err = zeros(cv.NumTestSets,1);
for i = 1:cv.NumTestSets
    trainIndex = allIndex(cv.training(i));

    % ignore indices that are not in the training set
    p = PARAMS;
    p.indexPos = PARAMS.indexPos(ismember(PARAMS.indexPos, trainIndex));
    p.indexNeg = PARAMS.indexNeg(ismember(PARAMS.indexNeg, trainIndex));
    
    % train HMM
    hmm = ForestHMM(data, p);
    
    fprintf(['\t> Testing samples...\n' ...
             '\t  ----------------------------\n']);
    % test on unseen data
    testIndex = allIndex(cv.test(i));
    
    ytest = zeros(numel(testIndex), 1);
    yclass = ismember(testIndex, PARAMS.indexPos); % the true class
    for j = 1:numel(testIndex)
        fprintf('Testing sample with index %d\n', testIndex(j));
        
        fprintf('SIZE');
        disp(size(data(:, :, testIndex(j))));
        ytest(j) = hmm.predict(data(:, :, testIndex(j)));
        fprintf('done');
    end
    err(i) = sum(~strcmp(ytest, yclass));
end
cvErr = sum(err)/sum(cv.TestSize);


% if PARAMS.method == 'forest'
%     hmm = ForestHMM(data, PARAMS);
%     hmm.predict(
% end




fprintf('Plotting results...\n');