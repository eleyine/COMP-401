function [ subsetStruct ] = getSubsetStruct( originalStruct, fieldSubset)
%GETSUBSETSTRUCT Returns struct with subset of fields

F = fieldnames(originalStruct);
D = struct2cell(originalStruct);
% make boolean mask with 1 corresponding to fields to keep
keep = zeros(numel(F));
for i = 1:numel(F)
    if find(strcmp(fieldSubset, F))
        keep(i) = 1;
    else
        keep(i) = 0;
    end
end

subsetStruct = cell2struct(D(keep), F(keep));
end
