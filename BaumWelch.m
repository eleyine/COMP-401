function [LL, prior1, transmat1, mu1, Sigma1, mixmat1] = BaumWelch(data, params)
%BAUMWELCH Estimation of parameters using Baum-Welch (EM)
    O = size(data, 1);          %Number of coefficients in a vector 
    T = size(data, 2);          %Number of vectors in a sequence 
    nex = size(data, 3);        %Number of sequences 
    M = params.M;               %Number of mixtures 
    Q = params.Q;               %Number of states 
    maxIter = params.maxIter;
    cov_type = 'full';

    % initial guess of parameters
    prior0 = normalise(rand(Q,1));
    transmat0 = mk_stochastic(rand(Q,Q));
    
    [mu0, Sigma0] = mixgauss_init(Q*M, data, cov_type);
    mu0 = reshape(mu0, [O Q M]);
    Sigma0 = reshape(Sigma0, [O O Q M]);
    mixmat0 = mk_stochastic(rand(Q,M));

    [LL, prior1, transmat1, mu1, Sigma1, mixmat1] = ...
    mhmm_em(data, prior0, transmat0, mu0, Sigma0, mixmat0, 'max_iter', maxIter);
end

