classdef ForestHMM
    %ForestHMM: trains a single HMM for every sample (baby) and classifies
    %an unknown sequence of observation based on the average log likelihood of
    %negative vs. positive classifiers
    
    properties
      pHMMs = struct;
      nHMMs = struct;
      indexPos = [];
      indexNeg = [];
    end
    
    methods
        function obj = ForestHMM(data, params)
            fprintf('*** BUILDING FOREST HMM ***\n');
            fprintf(['\t> Training Positive samples...\n' ...
                    '\t  ----------------------------\n']);
            obj.indexPos = params.indexPos;
            obj.indexNeg = params.indexNeg;
            p = numel(obj.indexPos);
            
            % preallocate p positive HMMS
            obj.pHMMs = repmat(struct(...
                    'LL', [], ...
                    'prior', [], ...
                    'transmat', [], ...
                    'mu', [], ...
                    'sigma', [], ...
                    'mixmat', []), p, 1);
            for i = 1:p
                fprintf('\t\t * %d\n', obj.indexPos(i));
                d = data(:,:,obj.indexPos(i));
                disp(obj.pHMMs);
                fprintf('SIZE');
                disp(size(d));
                [obj.pHMMs(i).LL, obj.pHMMs(i).prior, obj.pHMMs(i).transmat, ...
                    obj.pHMMs(i).mu, obj.pHMMs(i).sigma, obj.pHMMs(i).mixmat] = ...
                    BaumWelch(d, params);
                fprintf('\n');
            end
            fprintf(['\t> Training Negative samples...\n' ...
                     '\t  ----------------------------\n']);
            n = numel(obj.indexNeg);
            % preallocate n negative HMMS
            obj.nHMMs = repmat(struct(...
                    'LL', [], ...
                    'prior', [], ...
                    'transmat', [], ...
                    'mu', [], ...
                    'sigma', [], ...
                    'mixmat', []), n, 1);
            for i = 1:n
                fprintf('\t\t * %d\n', obj.indexNeg(i));
                d = data(:,:,obj.indexNeg(i));
                disp(obj.nHMMs);
                [obj.nHMMs(i).LL, obj.nHMMs(i).prior, obj.nHMMs(i).transmat, ...
                    obj.nHMMs(i).mu, obj.nHMMs(i).sigma, obj.nHMMs(i).mixmat] = ...
                    BaumWelch(d, params);
                fprintf('\n');
            end
        end
        
        function classValue = predict(obj, sequence)
%            fprintf('*** PREDICTING SEQUENCE USING FOREST HMM ***\n');
           ploglik = -Inf(numel(obj.indexPos));
           for i = 1:numel(obj.indexPos)
               ploglik(i) = mhmm_logprob(sequence, obj.pHMMs(i).prior, ...
                   obj.pHMMs(i).transmat, obj.pHMMs(i).mu, ...
                   obj.pHMMs(i).sigma, obj.pHMMs(i).mixmat);
           end
            nloglik = -Inf(numel(obj.indexNeg));
            for i = 1:numel(obj.indexNeg)
                nloglik(i) = mhmm_logprob(sequence, obj.nHMMs.prior(i), ...
                    obj.nHMMs.transmat(i), obj.nHMMs.mu(i), ...
                    obj.nHMMs.sigma(i), obj.nHMMs.mixmat(i));
            end
%            fprintf('p loglikelihood: %d vs n loglikelihood: %d\n', mean(ploglik), mean(nloglik));
            fprintf('==> %d\n', mean(ploglik) > mean(nloglik));
           classValue = mean(ploglik) > mean(nloglik);
        end
    end
    
end

