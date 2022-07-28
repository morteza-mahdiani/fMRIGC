classdef Connectivity
    %Summary of this class
    properties
        inputPathOfData
        inputPathOfMasks
        outputPath
        first_sub_ID
        last_sub_ID
%         prepDataList = containers.Map
    end
    properties
        sub_range
        regions
    end

    methods
        function obj = Connectivity(pathOfData, pathOfMasks, outputPath, fID, lID, ROIs, path_to_MVGC)
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
            obj.regions = ROIs;

            obj.setPath(path_to_MVGC);

        end
        function preprocess(obj, TCDataLength, saveFlag)
            seedStr = '%s_%s_10mm_Sphere.nii';

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

                allSeedTCMat = nan(TCDataLength,size(obj.regions,2));
                for sd =1:size(obj.regions,2)
                    cSubjSeed = fullfile(obj.inputPathOfMasks,cSubj, sprintf(seedStr,cSubj,obj.regions{sd}));
                    if ~exist(cSubjSeed),continue,end

                    % use cosmo to load data
                    ds_seed = cosmo_fmri_dataset(cSubjTC,'mask',cSubjSeed);
                    ds_seed = cosmo_remove_useless_data(ds_seed);
                    allSeedTCMat(:,sd) = mean(ds_seed.samples,2);
                end
                % store preprocessed data
%                 obj.prepDataList{end + 1} = transpose(allSeedTCMat);


                % save the preprocessed data if save flag is true
                if saveFlag == true
                    save(append(obj.outputPath, append(cSubj, '.mat')), 'obj.prepData');
                end
            end
        end

        function GCM(obj, GC_param_obj, actual_model_order, preprocessed_data_path)
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

            GC3DMat = nan(length(obj.regions), length(obj.regions), length(obj.sub_range));

            for n = 1 : length(listOfsubjs)
                %% Parameters
                ntrials   = GC_param_obj.ntrials;     % number of trials
                nobs      = GC_param_obj.nobs;   % number of observations per trial
                regmode   = GC_param_obj.regmode;  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
                icregmode = GC_param_obj.icregmode;  % information criteria regression mode ('OLS', 'LWR' or empty for default)
                morder    = GC_param_obj.morder;  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
                momax     = GC_param_obj.momax;     % maximum model order for model order estimation
                tstat     = GC_param_obj.tstat;     % statistical test for MVGC:  'chi2' for Geweke's chi2 test (default) or'F' for Granger's F-test
                alpha     = GC_param_obj.alpha;   % significance level for significance test
                mhtc      = GC_param_obj.mhtc;  % multiple hypothesis test correction (see routine 'significance')
                seed      = GC_param_obj.seed;      % random seed (0 for unseeded)

                %% Generate VAR test data
                nvars = length(obj.regions); % number of variables

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

        function GCTensor(obj, GC_param_obj, actual_model_order, preprocessed_data_path, loadFlag)
            if nargin == 2
                preprocessed_data_path = nan;
                loadFlag = false;
            elseif nargin == 3
                loadFlag = true;
            end

            if loadFlag == true
                obj.GCM(GC_param_obj, actual_model_order, preprocessed_data_path)
            elseif loadFlag == false
                %obj.GCM
            end

        end
        function visualize(obj, path_to_data, mode, subjectID)
            if isfile(path_to_data)
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
