function [ subsamples ] = subsample( features, params )
%SUBSAMPLE Subsample the features according to a params with fields
%   noNaN: boolean
%       Specifies whether or not to remove leading and trailing NaN from
%       feature columns.
%   timeFrame: [start end]
%       The timeframe in minutes from which to trim the arrays
%      '': 
% 
fields = fieldnames(features);
subsamples = struct(features);
fprintf('\t\tSample size now:\n'); disp(subsamples);

% Ignore some samples
if numel(params.ignoreSamples) > 0
    fprintf('\t> Removing some baby samples...\n');
    assert(max(params.ignoreSamples) <= numel(features.(fields{1})));
    keep = 1:numel(features.(fields{1}));
    keep(params.ignoreSamples) = [];
    for i = 1:numel(fields)
        subsamples.(fields{i}) = features.(fields{i})(:,keep);
    end
    fprintf('\t\t Sample size now:\n'); disp(subsamples);
end


% Remove trailing NaN from arrays
if params.noNaN == 1
    fprintf('\t> Removing Leading and Trailing NaN from feature columns...\n');
    for i = 1:numel(fields)
        f = subsamples.(fields{i});
        for j = 1:size(f, 2)
            sample = f{j};
            s = find(~isnan(sample), 1, 'first'); % start
            e = find(~isnan(sample), 1, 'last'); % end
            subsamples.(fields{i}){j} = sample(s:e);
            % Assert there are no more NaNs in the data
            assert(sum(isnan(subsamples.(fields{i}){j})) == 0);
        end
    end
    fprintf('\t\tSample size now:\n'); disp(subsamples);
end

% Reformatting structure of dimension (features, samples, timeframes) to a
% matrix of dimension (features, timeframes, samples) that is more
% convenient for HMM analysis
fprintf('\t> Reformating data...\n');
subsamples = reformatData(subsamples);
fprintf('\t\tSample size now:'); disp(size(subsamples));

if numel(params.timeFrame) > 0
    % Trimming feature columns to equal lengths using specified timeframe
    fprintf('\t> Subsampling according to specified timeframe in minutes...\n');
    % Convert minutes to timeframe index considering the sampling frequency
    s = params.timeFrame(1) * 60 * params.sf;
    e = params.timeFrame(2) * 60 * params.sf;
    if s >=0 && e >=0
        subsamples = subsamples(:,s:e,:);
    elseif s <0 && e >=0
        subsamples = subsamples(:,end+s+1:e);
    elseif s >=0 && e <0
        subsamples = subsamples(:,s:end+e+1,:);
    else
        subsamples = subsamples(:,end+s+1:end+e+1,:);
    end
    fprintf('\t\tSample size now:'); disp(size(subsamples));
end


if params.samplingTime > 0
    % Subsampling every kth second
    fprintf('\t> Subsampling every kth second as specified...\n');
    % find corresponding index
    k = floor(params.samplingTime * params.sf);
    subsamples = subsamples(:, 1:k:end,:);
    fprintf('\t\tSample size now:'); disp(size(subsamples));
end

% Fill remaining NaNs
s = sum(sum(sum(isnan(subsamples))));
if s > 0
    fprintf('WARNING: your samples contain %i NaNs. Interpolating...', s);
    for i = 1:size(subsamples,3)
       if sum(sum(isnan(subsamples(:,:,i))))> 0
         subsamples(:,:,i) = inpaint_nans(subsamples(:,:,i));
       end
    end
    sum(sum(sum(isnan(subsamples))))
end

end
