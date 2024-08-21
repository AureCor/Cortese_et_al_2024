% This script calculates weights for iterative SLR for fMRI data
% Run job_preprocess_mrVista to obtain preprocessed dataset for each ROI.

clearvars;

% specify project directory
patternDir = './data_preproc/';
inductionDir = './results/';

% sbj information
sbjs = {'Sub03' 'Sub05' 'Sub06' 'Sub07' 'Sub08' 'Sub09' ...
    'Sub10' 'Sub11' 'Sub12' 'Sub13' 'Sub14' 'Sub15' 'Sub16' 'Sub17' 'Sub18' ...
    'Sub19' 'Sub20' 'Sub21' 'Sub22' 'Sub23' 'Sub24' 'Sub25' 'Sub26' 'Sub27' ...
    'Sub28' 'Sub29' 'Sub30' 'Sub31' 'Sub32' 'Sub33' 'Sub34' 'Sub35' 'Sub36' ...
    'Sub37' 'Sub38' 'Sub39' 'Sub40' 'Sub41' 'Sub42' 'Sub43' 'Sub44'};

Nsbj    = length(sbjs);

rois    = {'ACC_pre_R', 'ACC_sub_L', 'ACC_sub_R', 'ACC_sup_L', 'ACC_sup_R', 'Amygdala_L', 'Amygdala_R', 'Angular_L', 'Angular_R', 'Calcarine_L', 'Calcarine_R',...
    'Caudate_L','Caudate_R', 'Cerebellum_10_L', 'Cerebellum_10_R', 'Cerebellum_3_L', 'Cerebellum_3_R', 'Cerebellum_4_5_L', 'Cerebellum_4_5_R', 'Cerebellum_6_L',...
    'Cerebellum_6_R', 'Cerebellum_7b_L', 'Cerebellum_7b_R', 'Cerebellum_8_L', 'Cerebellum_8_R', 'Cerebellum_9_L', 'Cerebellum_9_R', 'Cerebellum_Crus1_L', 'Cerebellum_Crus1_R',...
    'Cerebellum_Crus2_L', 'Cerebellum_Crus2_R', 'Cingulate_Mid_L', 'Cingulate_Mid_R', 'Cingulate_Post_L', 'Cingulate_Post_R', 'Cuneus_L', 'Cuneus_R', 'Frontal_Inf_Oper_L',...
    'Frontal_Inf_Oper_R', 'Frontal_Inf_Orb_2_L', 'Frontal_Inf_Orb_2_R', 'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Med_Orb_L', 'Frontal_Med_Orb_R', 'Frontal_Mid_2_L',...
    'Frontal_Mid_2_R', 'Frontal_Sup_2_L', 'Frontal_Sup_2_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'Fusiform_L', 'Fusiform_R', 'Heschl_L', 'Heschl_R', 'Hippocampus_L',...
    'Hippocampus_R', 'Insula_L', 'Insula_R', 'Lingual_L', 'Lingual_R', 'N_Acc_L', 'N_Acc_R', 'OFCant_L', 'OFCant_R', 'OFClat_L', 'OFClat_R','OFCmed_L', 'OFCmed_R', 'OFCpost_L',...
    'OFCpost_R', 'Occipital_Inf_L', 'Occipital_Inf_R', 'Occipital_Mid_L', 'Occipital_Mid_R', 'Occipital_Sup_L', 'Occipital_Sup_R', 'Olfactory_L', 'Olfactory_R', 'Pallidum_L',...
    'Pallidum_R', 'ParaHippocampal_L', 'ParaHippocampal_R', 'Paracentral_Lobule_L', 'Paracentral_Lobule_R', 'Parietal_Inf_L', 'Parietal_Inf_R', 'Parietal_Sup_L', 'Parietal_Sup_R',...
    'Postcentral_L', 'Postcentral_R', 'Precentral_L', 'Precentral_R', 'Precuneus_L', 'Precuneus_R', 'Putamen_L', 'Putamen_R', 'Rectus_L', 'Rectus_R', 'Rolandic_Oper_L',...
    'Rolandic_Oper_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R', 'SupraMarginal_L', 'SupraMarginal_R', 'Temporal_Inf_L', 'Temporal_Inf_R', 'Temporal_Mid_L', 'Temporal_Mid_R',...
    'Temporal_Pole_Mid_L', 'Temporal_Pole_Mid_R', 'Temporal_Pole_Sup_L', 'Temporal_Pole_Sup_R', 'Temporal_Sup_L', 'Temporal_Sup_R', 'Vermis_10', 'Vermis_1_2', 'Vermis_3',...
    'Vermis_4_5', 'Vermis_6', 'Vermis_7', 'Vermis_8', 'Vermis_9', 'Thal_L', 'Thal_R', 'VTA', 'SN', 'Raphe', 'Red', 'LC'};


% class label.
classLabel = [1 2];
% specify the number of iteration for SLR learning.
nLearn = 50;
% specify the maximum number of SLR.
Nslr = 10;
Ncv  = 20;
trte = 0.8; % train - test split

% parameters
driftRejection  = 1;
baselinetype    = 3; %1: Pre-blank, 2: Current data
minsample       = 8;

BOLDdelay       = 3;

session = {'Con' 'Ext' 'Tes'}; % can be: 'Con' (conditioning), 'Ext' (extinction), 'Test' (test)
cstype  = [1 2]; % S-cs+ vs E-cs+
cslabel = {'scspecsp'};

%%%%%%%%%%%%%%%%%%%%%%%
% process starts here %
%%%%%%%%%%%%%%%%%%%%%%%
for ses=1:length(session)
    for sub=1:Nsbj
        for ROI_num=1:length(rois) % repeate the following process for each ROI
            rng(3213); % set for reproducibility
            
            % sub
            Sb = sbjs{sub};
            
            
            clear *data* Ptr* Pte* dummy*;
            % Following cmb for detection task data
            % specify file postfix according to preprocess parameters
            postfix = ['_' session{ses} '_fMRI_lin_zscore_wbaseln_bd' num2str(BOLDdelay)];
            
            roiName = rois{ROI_num};
            disp(['Sub: ' Sb])
            disp(['Spatiotemporal iSLR for ' roiName ' ++++++++++++++++']);
            
            % load ROI data
            if exist([patternDir Sb '/rin_' roiName postfix '.mat'],'file')
                roidat = load([patternDir Sb '/rin_' roiName postfix '.mat'],'data');
                
                
                %%%%%%%%%%%%%%%%%%%
                % Change here to select the data set
                datatrain   = roidat.data.fullsoundseq(cstype)';
                datatest    = [roidat.data.seq3'; roidat.data.us3'];
                dect        = 'soundvsseq';
                Nclass      = size(datatrain,1);
                tmp         = cell2mat(cellfun(@size, datatrain, 'UniformOutput', false));
                Ntrialtr    = tmp(:,2);
                Nvoxel      = unique(tmp(:,1));
                %%%%%%%%%%%%%%%%%%%%
                
                
                % start process only if data samples are large enough for each
                % class
                if sum(cellfun('size',datatrain,2)>minsample)==Nclass
                    
                    % error if class number if wrong
                    if size(datatrain,1)~=length(classLabel)
                        error('Invalid class number');
                    end
                    
                    
                    % prepare train, test data and labels with data from sound
                    % sequence. This will be used to compute the optimal number
                    % of SLR. Then, we train the decoder with all data from the
                    % sound sequences, and test it on the data from the truck
                    % scene with optimal number of SLR runs.
                    X = cell2mat(datatrain')'; label = [];
                    for c = 1:Nclass
                        label = [label; c*ones(Ntrialtr(c),1)];
                    end
                    
                    
                    % initialize variables based on data dimensions. The first
                    % dimension of weights should be Nvoxel + 1 to inlude the bias term
                    % calculated by SMLR
                    weights_islr   = zeros(Nvoxel+1, Ncv, Nslr);
                    acctrm_islr    = [];
                    acctem_islr    = [];
                    acctr_islr     = [];
                    accte_islr     = [];
                    Nfeatsel_islr  = [];
                    Nvoxsel_islr   = [];
                    trtsample       = sum(floor(Ntrialtr*trte));
                    tetsample       = sum(Ntrialtr) - trtsample;
                    trsample        = floor(Ntrialtr*trte);
                    tesample        = Ntrialtr - trsample;
                    Ptr             = zeros(trtsample, Nclass, Nslr, Ncv);
                    Pte             = zeros(tetsample, Nclass, Nslr, Ncv);
                    Ptrm_islr      = zeros(trtsample, Nclass, Nslr, Ncv);
                    Ptem_islr      = zeros(tetsample, Nclass, Nslr, Ncv);
                    
                    % repeat cross-validation Ncv times
                    for CV=1:Ncv
                        disp(['-- Leave 1 out cross-validation to determine optimal number of SLRs: ' sprintf('%d',CV)]);
                        
                        [ixtr,ixte] = separate_train_test(label, trte);
                        % assign labels
                        ttr = label(ixtr,:);
                        tte = label(ixte,:);
                        
                        % keep track of voxels that have a weight assigned
                        nonzeros = [];
                        
                        % repeat Nslr times to find the optimal number of
                        % iterations
                        for SLR=1:Nslr
                            disp(['------ #SLR: ' sprintf('%d', SLR)]);
                            index = 1:Nvoxel;
                            index(nonzeros) = []; % remove voxels for which decoder weight has been assigned
                            
                            disp(['--------  ' sprintf('%d', length(index)) ' voxels for training']);
                            
                            
                            % use voxels which do not have decoder weight yet
                            % first, select data for train/test according to
                            % this CV run (same as labels defined outside this
                            % loop)
                            xtr = X(ixtr,:);
                            xte = X(ixte,:);
                            % second, select only voxels for which decoder
                            % weight has yet to be assigned
                            xtr = xtr(:,index);
                            xte = xte(:,index);
                            
                            
                            % Binary classification with SLR
                            [ww, ix_eff, errTable_tr, errTable_te, parm, AXall, Ptr(:, :, SLR, CV), Pte(:, :, SLR, CV)] = ...
                                biclsfy_slrvar_bal(xtr, ttr, xte, tte, 'nlearn', nLearn, 'mean_mode', 'none', 'scale_mode', 'none');
                            
                            % training accuracy for this SLR
                            acctr_islr(CV,SLR) = mean(diag(errTable_tr)./trsample);
                            % test accuracy for this SLR
                            accte_islr(CV,SLR) = mean(diag(errTable_te)./tesample);
                            % weight for this SLR
                            weights_islr([index Nvoxel+1],CV,SLR) = ww;
                            % number of selected features (basically, the sum
                            % of selected voxels across classes - can contain
                            % voxels selected twice)
                            Nfeatsel_islr(CV,SLR) = length(find(abs(reshape(ww(1:length(index),:),length(index)*1,1))>0));
                            % number of selected voxels (unique voxels)
                            Nvoxsel_islr(CV,SLR) = length(find(any(ww,2)==1));
                            % find voxels for which weight has been assigned in
                            % this SLR calculation
                            nonzeros = [nonzeros; find(any(weights_islr(1:Nvoxel,CV,SLR),2)==1)];
                            % error table (confusion matrix)
                            errorTable(:,:,SLR,CV) = errTable_te;
                            
                            
                            % calculate class probability for each class and SLR
                            if SLR==1
                                dummy1(:,:,SLR) = Ptr(:,:,SLR,CV);
                                dummy2(:,:,SLR) = Pte(:,:,SLR,CV);
                            else
                                % output of iSLR is calculated as multiplication from
                                dummy1(:,:,SLR) = dummy1(:,:,SLR-1) .* Ptr(:,:,SLR,CV);
                                dummy2(:,:,SLR) = dummy2(:,:,SLR-1) .* Pte(:,:,SLR,CV);
                            end
                            % normalize class probabilities
                            for CLASS=1:Nclass
                                Ptrm_islr(:,CLASS,SLR,CV) = dummy1(:,CLASS,SLR) ./ sum(dummy1(:,:,SLR),2);
                                Ptem_islr(:,CLASS,SLR,CV) = dummy2(:,CLASS,SLR) ./ sum(dummy2(:,:,SLR),2);
                            end
                            % calculate training and test accuracies of iSLR for each class
                            for CLASS=1:Nclass
                                classIndex = find(ttr==CLASS);
                                [~,imax] = max(Ptrm_islr(classIndex,:,SLR,CV),[],2);
                                acctrm_islr(CV,SLR,CLASS) = sum(imax==classLabel(CLASS))/length(imax);
                                
                                classIndex = find(tte==CLASS);
                                [pmax,imax] = max(Ptem_islr(classIndex,:,SLR,CV),[],2);
                                acctem_islr(CV,SLR,CLASS) = sum(imax==classLabel(CLASS))/length(imax);
                            end
                            
                            % display results for this SLR
                            disp('  ')
                            disp(['--------  ' sprintf('%d', Nfeatsel_islr(CV,SLR)) ' features selected']);
                            disp(['--------  ' sprintf('%d', Nvoxsel_islr(CV,SLR)) ' voxels selected']);
                            disp(['--------  training accuracy: ' sprintf('%d', round(acctr_islr(CV,SLR)*100)) '%']);
                            disp(['--------  test accuracy: ' sprintf('%d', round(accte_islr(CV,SLR)*100)) '%']);
                            disp('  ')
                        end
                    end
                end
                
                % save all results
                if ~exist([inductionDir Sb],'dir')
                    mkdir([inductionDir Sb]);
                end
                save([inductionDir Sb '/rin_' roiName postfix '_' dect '_' cslabel{cstype(1)} '_ibislr_weights_train.mat'], 'weights*', 'ww*', 'acc*','N*','err*', 'nonzeros');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                % Use the SLR yielding highest classification accuracy
                [accmax, Nslropt] = max(squeeze(mean(mean(acctem_islr,3),1)));
                
                disp(['------ optimal #SLR: ' sprintf('%d', Nslropt)]);
                disp(['------ test accuracy for optimal #SLR: ' sprintf('%d', round(accmax*100)) '%'])
                
                % prepare training data
                Xtr = cell2mat(datatrain')'; ttr = [];
                for c = 1:Nclass
                    ttr = [ttr; c*ones(Ntrialtr(c),1)];
                end
                
                % prepare test data
                tmp         = cell2mat(cellfun(@size, datatest, 'UniformOutput', false));
                Ntrialte    = tmp(cstype,2);
                Xte = cell2mat(datatest(cstype)')'; tte = [];
                for c = 1:Nclass
                    tte = [tte; c*ones(Ntrialte(c),1)];
                end
                
                % initialize empty matrices and clear some old variables
                ww   = zeros(Nvoxel+1, Nslropt);
                Pte1 = []; clear dummy* acc* err*
                nonzeros = [];
                
                
                for SLR=1:Nslropt
                    disp(['------ #SLR: ' sprintf('%d', SLR)]);
                    index = 1:Nvoxel;
                    index(nonzeros) = []; % remove voxels for which decoder weight has been assigned
                    
                    disp(['--------  ' sprintf('%d', length(index)) ' voxels for training']);
                    
                    
                    
                    % use only voxels for which decoder weight has yet to be
                    % assigned
                    xtr = Xtr(:,index);
                    xte = Xte(:,index);
                    
                    [weights, ix_eff1, errTable_tr1, errTable_te1, ~, ~, ~, Pte1(:,:,SLR)] = ...
                        biclsfy_slrvar_bal(xtr, ttr, xte, tte, 'nlearn', nLearn, 'mean_mode', 'none', 'scale_mode', 'none');
                    
                    % weight for this SLR
                    ww([index Nvoxel+1],SLR) = weights; % the same weights are selected in both cases
                    % find voxels for which weight has been assigned in
                    % this SLR calculation
                    nonzeros = [nonzeros; find(any(ww(1:Nvoxel,SLR),2)==1)];
                    
                    
                    
                    % calculate class probability for each class and SLR
                    if SLR==1
                        dummy1(:,:,SLR) = Pte1(:,:,SLR);
                    else
                        % output of iSLR is calculated as multiplication from
                        dummy1(:,:,SLR) = dummy1(:,:,SLR-1) .* Pte1(:,:,SLR);
                    end
                    % normalize class probabilities
                    for CLASS=1:Nclass
                        Pte1m(:,CLASS,SLR) = dummy1(:,CLASS,SLR) ./ sum(dummy1(:,:,SLR),2);
                    end
                    % calculate training and test accuracies of iSLR for each class
                    for CLASS=1:Nclass
                        classIndex = find(tte==CLASS);
                        [~,imax] = max(Pte1m(classIndex,:,SLR),[],2);
                        accte1(SLR,CLASS) = sum(imax==classLabel(CLASS))/length(imax);
                    end
                end
                
                % calculate accuracy at last iteration
                acc_gen = diag(errTable_te1)./Ntrialte(1:2,1);
                
                
                % select voxels with non-zero weight
                w_tmp = weights_islr(1:end-1,:,1:Nslropt);
                if Nslropt>1
                    w_tmp = reshape(w_tmp, [size(w_tmp,1) size(w_tmp,2)*size(w_tmp,3)]);
                end
                
                clear nonzero_index weights
                nonzero_index = find(sum(abs(w_tmp')~=0));
                weights = w_tmp(sort(nonzero_index),:);
                
                % update postfix (data from 3TRs --> truck scene + ITI)
                postfix = [postfix '_ts3'];
                
                % save pruned weights file for each ROI
                save([inductionDir Sb '/rin_' roiName postfix '_' dect '_' cslabel{cstype(1)} '_ibislr_weights_test.mat'], 'weights', 'ww', 'nonzero*', 'acc*', 'err*');
                
            end
        end
    end
end
