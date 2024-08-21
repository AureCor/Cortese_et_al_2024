% job_SLiR
%
% This script runs sparse linear regression analysis of fMRI data.
% Regressor performance is evaluated through cross-validation.
%
% Aurelio Cortese 2022/08/10

clearvars;

% directory information
patternDir  = './data_preproc/';
resultDir   = './results/';

% subject info
subjects = {'Sub03' 'Sub05' 'Sub06' 'Sub07' 'Sub08' 'Sub09' ...
    'Sub10' 'Sub11' 'Sub12' 'Sub13' 'Sub14' 'Sub15' 'Sub16' 'Sub17' 'Sub18' ...
    'Sub19' 'Sub20' 'Sub21' 'Sub22' 'Sub23' 'Sub24' 'Sub25' 'Sub26' 'Sub27' ...
    'Sub28' 'Sub29' 'Sub30' 'Sub31' 'Sub32' 'Sub33' 'Sub34' 'Sub35' 'Sub36' ...
    'Sub37' 'Sub38' 'Sub39' 'Sub40' 'Sub41' 'Sub42' 'Sub43' 'Sub44'};

Nsbj    = length(subjects);

sessions = {'Con' 'Ext' 'Tes'};

% SLiR parameter
parm.Ntrain = 2000; % # of total training iteration
parm.Nskip = 100;	% skip steps for display info
parm.Tau = 0;   % Lag time steps
parm.Dtau = 1;   % Number of embedding dimension
parm.Tpred = 0;   % Prediction time step : y(t+Tpred) = W * x(t)
parm.data_norm = 1;	% Normalize input and output

% scaling parameter
linearize = 1; % 0 for no linearization, 1 for likelihood linearization using atanh
correction = 10e-16;

% roi info
sROI = {'HPC', 'BA46_9'}; % source ROI (likelihood)
tROI = {'vmPFC', 'amyg', 'amyg-vmPFC'}; % target ROI (timecourse)
% tROI = {'Insula', 'thal', 'fusiform', 'auditorycortex'}; % target ROI (timecourse)

Nclass = 2; % S-cs+ & E-cs+

disp('================');
disp('SLiR started');
disp('================');
for SBJ=1:length(subjects)
    
    disp(['Subject: ' subjects{SBJ}]);
    
    for sROI_index=1:length(sROI)
        for tROI_index=1:length(tROI)
            sRd = load([patternDir subjects{SBJ} '/rin_' sROI{sROI_index} '_trialwise_likel_roidata.mat'],'likelihood','trials');
            tRd = load([patternDir subjects{SBJ} '/rin_' tROI{tROI_index} '_trialwise_likel_roidata.mat'],'roidata');
            disp([sROI{sROI_index} ' to ' tROI{tROI_index}]);
            
            clear data result SLiRModel parameter
            
            for exps=1:length(sessions)
                
                for cc=1:Nclass
                    
                    Nrun = sRd.trials{exps}(1); % symmetric number of trials
                    
                    for run=1:Nrun
                        disp(['-- CV run ' sprintf('%d', run)]);
                        
                        % select training and test data
                        trainrun = 1:Nrun;
                        trainrun(run) = [];
                        testtrial = [run Nrun+run];
                        traintrial = 1:Nrun*2;
                        traintrial(testtrial) = [];
                        
                        xtrain = tRd.roidata{exps}(:,traintrial);
                        ytrain = sRd.likelihood{exps}(cc,traintrial);
                        
                        xtest = tRd.roidata{exps}(:,testtrial);
                        ytest = sRd.likelihood{exps}(cc,testtrial);
                        
                        if linearize==1
                            ytrain = (ytrain+correction)/(1+2*correction);
                            ytrain = atanh(2*(ytrain-0.5));
                            ytest = (ytest+correction)/(1+2*correction);
                            ytest = atanh(2*(ytest-0.5));
                        end
                        
                        [Nvoxel Ntrain] = size(xtrain);
                        xtrain = reshape(xtrain,Nvoxel,1,Ntrain);
                        ytrain = reshape(ytrain,1,1,Ntrain);
                        
                        [Nvoxel Ntest] = size(xtest);
                        xtest = reshape(xtest,Nvoxel,1,Ntest);
                        ytest = reshape(ytest,1,1,Ntest);
                        
                        % normalize
                        [xtrain,nparm] = normalize_data(xtrain,parm.data_norm);
                        parm.xmean = nparm.xmean;
                        parm.xnorm = nparm.xnorm;
                        xtest = normalize_data(xtest,parm.data_norm,parm);
                        
                        [ytrain,nparm] = normalize_data(ytrain,parm.data_norm);
                        parm.ymean = nparm.xmean;
                        parm.ynorm = nparm.xnorm;
                        
                        % run SLiR
                        Model = [];
                        [Model,Info] = linear_sparse_space(xtrain,ytrain,Model,parm);
                        
                        % evaluation
                        ypred = predict_output(xtrain,Model,parm);
                        trainerr = sum((ytrain(:)-ypred(:)).^2)/sum(ytrain(:).^2);
                        dummycorr = corrcoef(ypred(:),ytrain(:));
                        traincorr = dummycorr(1,2);
                        data{exps}{run,cc}.pred_train = ypred;
                        
                        ypred = predict_output(xtest,Model,parm);
                        testerr = sum((ytest(:)-ypred(:)).^2)/sum(ytest(:).^2);
                        dummycorr = corrcoef(ypred(:),ytest(:));
                        testcorr = dummycorr(1,2);
                        data{exps}{run,cc}.pred_test = ypred;
                        
                        SLiRModel{exps}{run,cc} = Model;
                        data{exps}{run,cc}.xtrain = xtrain;
                        data{exps}{run,cc}.ytrain = ytrain;
                        result.trainerr{exps}(run,cc,:) = trainerr;
                        result.traincorr{exps}(run,cc,:) = traincorr;
                        data{exps}{run,cc}.xtest = xtest;
                        data{exps}{run,cc}.ytest = ytest;
                        result.testerr{exps}(run,cc,:) = testerr;
                        result.testcorr{exps}(run,cc,:) = testcorr;
                        parameter{exps}{run,cc} = parm;
                    end
                end
            end
            
            % save the data
            if linearize==1
                save([resultDir subjects{SBJ} '/rin_' sROI{sROI_index} '_to_' tROI{tROI_index} '_linear_SLiR.mat'],'data','result','SLiRModel','parameter');
            else
                save([resultDir subjects{SBJ} '/rin_' sROI{sROI_index} '_to_' tROI{tROI_index} '_SLiR.mat'],'data','result','SLiRModel','parameter');
            end
            
        end
    end
end
