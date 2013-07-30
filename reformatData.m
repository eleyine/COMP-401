function [ samples ] = reformatData(feat)
%REFORMATDATA formats structure of features to matrix of samples
% Feature structure with dimension (features, samples, timeframes) is transformed into a
% matrix with dimension (features, timeframes, samples)
featMat = struct;
fields = fieldnames(feat);
for i = 1:size(fields)
    cell = feat.(fields{i});
    maxSize = max(cellfun(@numel,cell));    % Get the maximum vector size
    fcn = @(x) vertcat(nan(maxSize-numel(x),1), x);  % Create an anonymous function
    rmat = cellfun(fcn,cell,'UniformOutput',false);  % Pad each cell with NaNs
    featMat.(fields{i}) = cell2mat(rmat);
end
featMat = struct2cell(featMat);
samples = permute(cat(3, featMat{:}), [3 1 2]);
end

