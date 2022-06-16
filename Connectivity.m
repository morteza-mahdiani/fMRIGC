classdef Connectivity
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    % JW test comment
    properties
        inputPathOfData
    end
    properties
        inputPathOfMasks
    end
    properties
        outputPath
    end

    methods
        function obj = Connectivity(pathOfData, pathOfMasks, outputPath)
            obj.inputPathOfData = pathOfData;
            obj.inputPathOfMasks = pathOfMasks;
            obj.outputPath = outputPath;
        end
        function preprocess(obj)
            subjRange = [3:30];
            seedStr = '%s_%s_10mm_Sphere.nii';
            %regions = {'lPMTG','lMFus','lIPL'};
            regions = {'rOFA','rFFA','rSTSF'};
            TCDataLength = 1100;
            % prepare the proposed path
            for sbj = subjRange
                cSubj = sprintf('sub-%1.2d',sbj);
                cSubjTC = sprintf(obj.inputPathOfData,cSubj);
                files = dir(cSubjTC);
                flag = [files.isdir];
                if isempty(flag),continue,end
                fileName = files(~flag).name;
                cSubjTC = append(cSubjTC, fileName);
                if ~exist(cSubjTC),continue,end

                allSeedTCMat = nan(TCDataLength,size(regions,2));
                for sd =1:size(regions,2)
                    cSubjSeed = fullfile(obj.inputPathOfMasks,cSubj, sprintf(seedStr,cSubj,regions{sd}));
                    if ~exist(cSubjSeed),continue,end

                    % use cosmo to load data
                    ds_seed = cosmo_fmri_dataset(cSubjTC,'mask',cSubjSeed);
                    ds_seed = cosmo_remove_useless_data(ds_seed);
                    allSeedTCMat(:,sd) = mean(ds_seed.samples,2);
                end
                % save data
                data = transpose(allSeedTCMat);
                save(append(append(obj.outputPath, cSubj), '.mat'),'data')
            end
        end

        function GCM(obj)
            subjRange = [3:30];
            seedStr = '%s_%s_10mm_Sphere.nii';
            listOfsubjs = dir(append(obj.inputPathOfData, '*.mat'));
            GC3DMat = nan(3, 3, 26);
            for n = 1 : size(listOfsubjs)
                %% Parameters
                ntrials   = 1;     % number of trials
                nobs      = 1100;   % number of observations per trial
                regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
                icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
                morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
                momax     = 20;     % maximum model order for model order estimation
                tstat     = 'F';     % statistical test for MVGC:  'chi2' for Geweke's chi2 test (default) or'F' for Granger's F-test
                alpha     = 0.05;   % significance level for significance test
                mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')
                seed      = 0;      % random seed (0 for unseeded)

                %% Generate VAR test data
                nvars = 3; % number of variables
                % Residuals covariance matrix.
                SIGT = eye(nvars);
                fprintf('\n ______________________________________________________________________');
                fprintf('\n');
                disp(append(listOfsubjs(n).folder, '/', listOfsubjs(n).name));
                fprintf('\n');
                sbj = append(listOfsubjs(n).folder, '/', listOfsubjs(n).name);
                ptic('\n*** var_to_tsdata... ');
                load(sbj);
                X = data;
                ptoc;
                name = split(listOfsubjs(n).name, '.');
                %dirTosbj =  '/home/proactionlab/Documents/projects/connectivity/ToolsArea_outputs/' + string(name(1));
                dirTosbj =  '/home/proactionlab/Documents/projects/connectivity/FaceArea_outputs/' + string(name(1));
                mkdir(dirTosbj);
                %% Model order estimation
                % Calculate information criteria up to max model order
                ptic('\n*** tsdata_to_infocrit\n');
                [AIC,BIC] = tsdata_to_infocrit(X,momax,icregmode);
                ptoc('*** tsdata_to_infocrit took ');
                [~,bmo_AIC] = min(AIC);
                [~,bmo_BIC] = min(BIC);
                % Plot information criteria.
                figure(1); clf;
                plot((1:momax)',[AIC BIC]);
                saveas(gcf,dirTosbj+'/'+string(name(1))+'_IC.png')
                legend('AIC','BIC');
                amo = 20; % actual model order
                fprintf('\nbest model order (AIC) = %d\n',bmo_AIC);
                fprintf('best model order (BIC) = %d\n',bmo_BIC);
                fprintf('actual model order     = %d\n',amo);
                % Select model order
                if  strcmpi(morder,'actual')
                    morder = amo;
                    fprintf('\nusing actual model order = %d\n',morder);
                elseif strcmpi(morder,'AIC')
                    morder = bmo_AIC;
                    fprintf('\nusing AIC best model order = %d\n',morder);
                elseif strcmpi(morder,'BIC')
                    morder = bmo_BIC;
                    fprintf('\nusing BIC best model order = %d\n',morder);
                else
                    fprintf('\nusing specified model order = %d\n',morder);
                end
                %% Granger causality estimation
                % Calculate time-domain pairwise-conditional causalities. Return VAR parameters
                % so we can check VAR.
                ptic('\n*** GCCA_tsdata_to_pwcgc... ');
                [F,A,SIG] = GCCA_tsdata_to_pwcgc(X,morder,regmode); % use same model order for reduced as for full regressions
                ptoc;
                % Check for failed (full) regression
                assert(~isbad(A),'VAR estimation failed');
                % Check for failed GC calculation
                assert(~isbad(F,false),'GC calculation failed');
                % Check VAR parameters (but don't bail out on error - GCCA mode is quite forgiving!)
                rho = var_specrad(A);
                fprintf('\nspectral radius = %f\n',rho);
                if rho >= 1,       fprintf(2,'WARNING: unstable VAR (unit root)\n'); end
                if ~isposdef(SIG), fprintf(2,'WARNING: residuals covariance matrix not positive-definite\n'); end
                % Significance test using theoretical null distribution, adjusting for multiple
                % hypotheses.
                pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat);
                sig  = significance(pval,alpha,mhtc);
                % Plot time-domain causal graph, p-values and significance.
                figure(2); clf;
                subplot(1,3,1);
                plot_pw(F);
                disp(pval);
                title('Pairwise-conditional GC');
                subplot(1,3,2);
                plot_pw(pval);
                title('p-values');
                subplot(1,3,3);
                plot_pw(sig);
                title(['Significant at p = ' num2str(alpha)])
                fprintf(2,'\nNOTE: no frequency-domain pairwise-conditional causality calculation in GCCA compatibility mode!\n');
                saveas(gcf,dirTosbj+'/'+string(name(1))+'.png')
                GC3DMat( :,: ,n)= pval;
                n = n + 1;
                clear bmo_BIC bmo_AIC X figure(1) figure(2);
            end
            %save(append(append(obj.outputPath, 'GC3Doutput_toolsAreas'), '.mat'),'GC3DMat')
            %save(append(append(obj.outputPath, 'GC3Doutput_faceAreas'), '.mat'),'GC3DMat')
        end
    end
end
