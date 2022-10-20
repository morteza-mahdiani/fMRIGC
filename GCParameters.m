% GCParameters is a class for fMRIGC tool parameters.
% 
% Author: Morteza Mahdiani, 2022, m72.morteza@gmail.com
% Last update: Morteza Mahdiani, 19/10/2022, m72.morteza@gmail.com
% ____________________________________________________________________

classdef GCParameters
    % class for the parameters that are needed through granger causality analysis
    properties
        ntrials   = 1;      % number of trials
        regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
        icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
        morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
        momax     = 20;     % maximum model order for model order estimation
        tstat     = 'F';    % statistical test for MVGC:  'chi2' for Geweke's chi2 test (default) or'F' for Granger's F-test
        alpha     = 0.05;   % significance level for significance test
        mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
        seed      = 0;      % random seed (0 for unseeded)
    end

    methods
        %class constructor
        function obj = GCParameters(ntrials, regmode, icregmode, morder, momax, tstat, alpha, mhtc, seed)
            obj.ntrials   = ntrials;
            obj.regmode   = regmode;
            obj.icregmode = icregmode;
            obj.morder    = morder;
            obj.momax     = momax;
            obj.tstat     = tstat;
            obj.alpha     = alpha;
            obj.mhtc      = mhtc;
            obj.seed      = seed;
        end

    end
end