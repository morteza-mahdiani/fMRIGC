% fMRIGC is a tool for granger causality analysis in fMRI data.
%
% fMRIGC investigates the connectivity in time courses of
% the fMRI data in both regio-wise and voxel-wise analysis modes.
% 
% Author: Morteza Mahdiani, 2022, m72.morteza@gmail.com
% Last update: Morteza Mahdiani, 19/10/2022, m72.morteza@gmail.com
% ____________________________________________________________________

classdef FMRIGC
    % the body of main class for granger causality
    properties
        inputPathOfData         % path to the dataset we want to work on
        inputPathOfMasks        % path to the masks we want to apply on datasets
        outputPath              % output directory that will be used for storing results and mid-level processed data
        first_sub_ID            % ID of the first subject in ascending order
        last_sub_ID             % ID of the last subject in descending order
        prepDataList            % list of pre-processed data for further analysis
        lengthOfTC              % number of time-courses to be analysed for granger causality
        regionWise = false;     % flag to select region-wise or voxel-wise data analysis
        observations            % number of observations per trial

    end

    properties
        sub_range               % range including the first and last subject ID
        regions                 % regions that will be considerd for granger causality through region-wise mode
    end

    methods

        function obj = FMRIGC(pathOfData, pathOfMasks, outputPath, fID, lID, number_of_observations, path_to_MVGC, ...
                regionWise, ROIs)
            if nargin == 7
                obj.regions = nan;
                obj.regionWise = false;
            elseif nargin == 8
                disp('Enter the ROIs');
                return
            elseif nargin == 9
                obj.regionWise = regionWise;
                obj.regions = ROIs;
            else
                disp('Check the number of input arguments');
                return

            end

            obj.inputPathOfData = pathOfData;
            obj.inputPathOfMasks = pathOfMasks;
            obj.outputPath = outputPath;

            % check whether output directory exists or not
            obj.outputPath = fullfile(obj.outputPath, '/');
            if ~exist(obj.outputPath, 'dir')
                mkdir(obj.outputPath);
            end
            obj.first_sub_ID = fID;
            obj.last_sub_ID = lID;
            obj.sub_range = (obj.first_sub_ID: obj.last_sub_ID);
            obj.observations = number_of_observations;

            obj.prepDataList = containers.Map;
            obj.lengthOfTC = containers.Map;

            obj.setPath(path_to_MVGC);

        end

        function regionWisePreprocess(obj, saveFlag)
            if obj.regionWise == false
                disp('Region-wise flag has set to be false. Create the instance again and set it true!')
                return
            end
            obj.lengthOfTC('length') = size(obj.regions,2);

            % prepare the proposed path
            for sbj = obj.sub_range
                cSubj = sprintf('sub-%1.2d',sbj);
                cSubjTC = fullfile(obj.inputPathOfData,cSubj,'/');
                files = dir(cSubjTC);
                flag = [files.isdir];
                if isempty(flag),continue,end
                fileName = files(~flag).name;
                cSubjTC = append(cSubjTC, fileName);
                if ~exist(cSubjTC),continue,end

                allSeedTCMat = nan(obj.observations, size(obj.regions, 2));
                for sd =1:size(obj.regions,2)
                    cSubjSeed = fullfile(obj.inputPathOfMasks, cSubj, obj.regions{sd});
                    if ~exist(cSubjSeed),continue,end

                    % use cosmo to load data
                    ds_seed = cosmo_fmri_dataset(cSubjTC,'mask',cSubjSeed);
                    ds_seed = cosmo_remove_useless_data(ds_seed);
                    allSeedTCMat(:,sd) = mean(ds_seed.samples,2);
                end
                % store preprocessed data
                if saveFlag == false
                    obj.prepDataList(cSubj) = transpose(allSeedTCMat);

                end

                % save the preprocessed data if save flag is true
                if saveFlag == true
                    data = transpose(allSeedTCMat);
                    save(append(obj.outputPath, append(cSubj, '.mat')), 'data');

                end
            end
        end

        function voxelWisePreprocess(obj, saveFlag)
            if obj.regionWise == true
                disp('Region-wise flag has set to be true. Create the instance again and set it false!')
                return
            end

            % prepare the proposed path
            for sbj = obj.sub_range
                cSubj = sprintf('sub-%1.2d',sbj);
                cSubjTC = fullfile(obj.inputPathOfData,cSubj,'/');
                files = dir(cSubjTC);
                flag = [files.isdir];
                if isempty(flag),continue,end
                fileName = files(~flag).name;
                cSubjTC = append(cSubjTC, fileName);
                if ~exist(cSubjTC),continue,end

                % use cosmo to load data
                ds_seed = cosmo_fmri_dataset(cSubjTC);
                ds_seed = cosmo_remove_useless_data(ds_seed);

                allSeedTCMat = nan(size(ds_seed.samples,1),1);
                allSeedTCMat(:) = mean(ds_seed.samples,2);
                obj.lengthOfTC('length') = size(ds_seed.samples,2);

                % store preprocessed data
                if saveFlag == false
                    obj.prepDataList(cSubj) = transpose(allSeedTCMat);
                end

                % save the preprocessed data if save flag is true
                if saveFlag == true
                    data = transpose(allSeedTCMat);
                    save(append(obj.outputPath, append(cSubj, '.mat')), 'data');
                end
            end
        end
        
        function GCMLoad(obj, GC_param_obj, actual_model_order, preprocessed_data_path)

            pathOfData = fullfile(preprocessed_data_path, '/');

            % check whether the path directory exists or not
            if ~isfolder(pathOfData)
                return;
            end

            listOfsubjs = dir(append(pathOfData, '*.mat'));

            % check whether there is preprocessed data or not
            if isempty(listOfsubjs)
                return;
            end

            GC3DMat = nan(obj.lengthOfTC('length'), obj.lengthOfTC('length'), length(obj.sub_range));

            for n = 1 : length(listOfsubjs)
                %% Parameters
                ntrials   = GC_param_obj.ntrials;       % number of trials
                nobs      = obj.observations;           % number of observations per trial
                regmode   = GC_param_obj.regmode;       % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
                icregmode = GC_param_obj.icregmode;     % information criteria regression mode ('OLS', 'LWR' or empty for default)
                morder    = GC_param_obj.morder;        % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
                momax     = GC_param_obj.momax;         % maximum model order for model order estimation
                tstat     = GC_param_obj.tstat;         % statistical test for MVGC:  'chi2' for Geweke's chi2 test (default) or'F' for Granger's F-test
                alpha     = GC_param_obj.alpha;         % significance level for significance test
                mhtc      = GC_param_obj.mhtc;          % multiple hypothesis test correction (see routine 'significance')
                seed      = GC_param_obj.seed;          % random seed (0 for unseeded)

                %% Generate VAR test data
                nvars = obj.lengthOfTC('length'); % number of variables

                % Residuals covariance matrix.
                SIGT = eye(nvars);
                fprintf('\n');
                disp(append(listOfsubjs(n).folder, '/', listOfsubjs(n).name));
                fprintf('\n');
                sbj = append(listOfsubjs(n).folder, '/', listOfsubjs(n).name);
                ptic('\n*** var_to_tsdata... ');
                load(sbj);
                X = data;
                ptoc;
                %name = split(listOfsubjs(n).name, '.');
                %dirTosbj = obj.outputPath + '/GC_plots_and_outputs/' + string(name(1));
                %mkdir(dirTosbj);

                %% Model order estimation
                % Calculate information criteria up to max model order
                ptic('\n*** tsdata_to_infocrit\n');
                [AIC,BIC] = tsdata_to_infocrit(X,momax,icregmode);
                ptoc('*** tsdata_to_infocrit took ');
                [~,bmo_AIC] = min(AIC);
                [~,bmo_BIC] = min(BIC);
                % Plot information criteria.
                %                 figure(1); clf;
                %                 plot((1:momax)',[AIC BIC]);
                %                 saveas(gcf,dirTosbj+'/'+string(name(1))+'_IC.png')
                %                 legend('AIC','BIC');
                amo = actual_model_order; % actual model order
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
                %                 figure(2); clf;
                %                 subplot(1,3,1);
                %                 plot_pw(F);
                %                 disp(pval);
                %                 title('Pairwise-conditional GC');
                %                 subplot(1,3,2);
                %                 plot_pw(pval);
                %                 title('p-values');
                %                 subplot(1,3,3);
                %                 plot_pw(sig);
                %                 title(['Significant at p = ' num2str(alpha)])
                %                 fprintf(2,'\nNOTE: no frequency-domain pairwise-conditional causality calculation in GCCA compatibility mode!\n');
                %                 saveas(gcf,dirTosbj+'/'+string(name(1))+'.png')

                %convert p-vals to z score
                zscore = norminv(1 - pval );
                GC3DMat( :,: ,n)= zscore;

                %Free up used memory for other loops
                clear bmo_BIC bmo_AIC X figure(1) figure(2);
            end

            % save data
            outPath = fullfile(obj.outputPath, 'GCMOutput', '/');

            % check whether the path exists or not
            if ~exist(outPath, 'dir')
                mkdir(outPath);
            end
            save(append(outPath, 'GC3DMat.mat'),'GC3DMat');
        end

        function GCM(obj, GC_param_obj, actual_model_order )

            % check whether there is preprocessed data or not
            if isempty(obj.prepDataList)
                return;
            end

            numberOfSubjects = length(obj.sub_range);
            GC3DMat = nan(obj.lengthOfTC('length'), obj.lengthOfTC('length'), numberOfSubjects);

            counter = 1;
            for k = keys(obj.prepDataList)
                %% Parameters
                ntrials   = GC_param_obj.ntrials;     % number of trials
                nobs      = obj.observations;   % number of observations per trial
                regmode   = GC_param_obj.regmode;  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
                icregmode = GC_param_obj.icregmode;  % information criteria regression mode ('OLS', 'LWR' or empty for default)
                morder    = GC_param_obj.morder;  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
                momax     = GC_param_obj.momax;     % maximum model order for model order estimation
                tstat     = GC_param_obj.tstat;     % statistical test for MVGC:  'chi2' for Geweke's chi2 test (default) or'F' for Granger's F-test
                alpha     = GC_param_obj.alpha;   % significance level for significance test
                mhtc      = GC_param_obj.mhtc;  % multiple hypothesis test correction (see routine 'significance')
                seed      = GC_param_obj.seed;      % random seed (0 for unseeded)

                %% Generate VAR test data
                nvars = obj.lengthOfTC('length'); % number of variables

                % Residuals covariance matrix.
                SIGT = eye(nvars);
                fprintf('\n');
                disp(k);
                fprintf('\n');
                data = obj.prepDataList(k{1});
                ptic('\n*** var_to_tsdata... ');
                %                 load(sbj);
                X = data;
                ptoc;
                %name = split(listOfsubjs(n).name, '.');
                %dirTosbj = obj.outputPath + '/GC_plots_and_outputs/' + string(name(1));
                %mkdir(dirTosbj);

                %% Model order estimation
                % Calculate information criteria up to max model order
                ptic('\n*** tsdata_to_infocrit\n');
                [AIC,BIC] = tsdata_to_infocrit(X,momax,icregmode);
                ptoc('*** tsdata_to_infocrit took ');
                [~,bmo_AIC] = min(AIC);
                [~,bmo_BIC] = min(BIC);
                % Plot information criteria.
                %                 figure(1); clf;
                %                 plot((1:momax)',[AIC BIC]);
                %                 saveas(gcf,dirTosbj+'/'+string(name(1))+'_IC.png')
                %                 legend('AIC','BIC');
                amo = actual_model_order; % actual model order
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
                %                 figure(2); clf;
                %                 subplot(1,3,1);
                %                 plot_pw(F);
                %                 disp(pval);
                %                 title('Pairwise-conditional GC');
                %                 subplot(1,3,2);
                %                 plot_pw(pval);
                %                 title('p-values');
                %                 subplot(1,3,3);
                %                 plot_pw(sig);
                %                 title(['Significant at p = ' num2str(alpha)])
                %                 fprintf(2,'\nNOTE: no frequency-domain pairwise-conditional causality calculation in GCCA compatibility mode!\n');
                %                 saveas(gcf,dirTosbj+'/'+string(name(1))+'.png')

                %convert p-vals to z score
                zscore = norminv(1 - pval );

                GC3DMat( :,: ,counter)= zscore;
                counter = counter + 1;
                %                 if counter > numberOfSubjects
                %                     disp('Subject range is not consistent with number of subjects!')
                %                     return
                %                 end

                %Free up used memory for other loops
                clear bmo_BIC bmo_AIC X figure(1) figure(2);
            end

            % save data
            outPath = fullfile(obj.outputPath, 'GCMOutput', '/');

            % check whether the path exists or not
            if ~exist(outPath, 'dir')
                mkdir(outPath);
            end
            save(append(outPath, 'GC3DMat.mat'),'GC3DMat');
        end

        function GCTensor(obj, GC_param_obj, actual_model_order, preprocessed_data_path)
            if nargin == 3
                preprocessed_data_path = nan;
                loadFlag = false;
            elseif nargin == 4
                loadFlag = true;
            else
                disp('Check the number of input arguments');
                return
            end

            if loadFlag == true
                obj.GCMLoad(GC_param_obj, actual_model_order, preprocessed_data_path)
            elseif loadFlag == false
                obj.GCM(GC_param_obj, actual_model_order)
            end

        end

        function regionWiseVisualization(obj, path_to_data, mode, subjectID)
            if isfile(path_to_data) && obj.regionWise == 1
                % load data
                conData = load(path_to_data);
                dataForVisualization = conData.GC3DMat;

                % eliminate NaN and INF datapoints
                dataForVisualization(isinf(dataForVisualization) ...
                    | isnan(dataForVisualization))= 0;

                if mode == 'm'
                    circularGraph(mean(dataForVisualization, 3), 'Label',obj.regions);

                elseif mode == 'sn'
                    if (1 <= subjectID) && (subjectID <= size(dataForVisualization, 3))
                        circularGraph(dataForVisualization(:,:,subjectID), 'Label',obj.regions);
                    else
                        disp('Enter a valid subject ID!');
                        return
                    end
                else
                    disp("Mode should be 'm' or 'sn'!")
                    return
                end

            else
                disp("The intended file for visualizayion doesn't exist!");
                return
            end

        end

    end


    methods(Static)
        function setPath(path_to_MVGC)
            pathOfMVGC = fullfile(path_to_MVGC, '/');

            % Add mvgc root directory and appropriate subdirectories to path
            mvgc_root = fileparts(mfilename(pathOfMVGC)); % directory containing this file

            % essentials
            addpath(mvgc_root);
            addpath(fullfile(mvgc_root,'core'));
            addpath(fullfile(mvgc_root,'gc'));
            addpath(fullfile(mvgc_root,'gc','GCCA_compat'));
            addpath(fullfile(mvgc_root,'gc','subsample'));
            addpath(fullfile(mvgc_root,'stats'));
            addpath(fullfile(mvgc_root,'utils'));
            if ~fexists(@rng) || ~fexists(@randi) % legacy hack
                addpath(fullfile(mvgc_root,'utils','legacy'));
                if ~fexists(@rng),   addpath(fullfile(mvgc_root,'utils','legacy','rng'));   end
                if ~fexists(@randi), addpath(fullfile(mvgc_root,'utils','legacy','randi')); end
            end
            addpath(fullfile(mvgc_root,'demo'));
            addpath(fullfile(mvgc_root,'mex'));
            addpath(fullfile(mvgc_root,'experimental'));
            addpath(fullfile(mvgc_root,'docs')); % don't add the 'html' subdirectory

            % comment out for release
            % addpath(fullfile(mvgc_root,'testing'));
            % addpath(fullfile(mvgc_root,'maintainer'));

            fprintf('[mvgc startup] Added MVGC root directory %s and subdirectories to path\n',mvgc_root);
        end
    end

end
