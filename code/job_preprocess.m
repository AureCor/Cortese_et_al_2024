% This preprocessing script is designed for the data obtained from the
% experimental sessions of the study: high-order sequences in threat
% learning.
% The basic design of the task included a movie with 3 sound elements
% presented in different order. Only one sequence was predictive of a truck
% crash. Decoding is done with trials labeled according to sequence type,
% event type (first or last) and data from the final scene when the crash may happen.
% In short, this script extract activation patterns from predefined ROIs.
% Such activation patterns correspond to time points of interest in the
% trial sequence (sounds, truck scene, etc).

% Aurelio, first developed 24-01-2020, last modified 03-08-2024


clearvars;

% basic parameters
sbjs = {'Sub03' 'Sub05' 'Sub06' 'Sub07' 'Sub08' 'Sub09' ...
    'Sub10' 'Sub11' 'Sub12' 'Sub13' 'Sub14' 'Sub15' 'Sub16' 'Sub17' 'Sub18' ...
    'Sub19' 'Sub20' 'Sub21' 'Sub22' 'Sub23' 'Sub24' 'Sub25' 'Sub26' 'Sub27' ...
    'Sub28' 'Sub29' 'Sub30' 'Sub31' 'Sub32' 'Sub33' 'Sub34' 'Sub35' 'Sub36' ...
    'Sub37' 'Sub38' 'Sub39' 'Sub40' 'Sub41' 'Sub42' 'Sub43' 'Sub44'};

Nsbj    = length(sbjs);

% specify project directory
dataDir     = '../data/';
patternDir  = './data_preproc/';
condDir     = '../data/Behav/';
anatDir     = '../data/Anatomy/';
filepfx     = 'sub*ext*run*'; % options: conditioning, extinction, test
condpfx     = 'Extinction'; % options: Conditioning, Extinction, Test, Rest
idx         = [6 7]; % indices change for different types of conditions
% acquisition: [1 2 3 4]
% extinction:   [6 7]
% test:         [10 11]
% rest:         [5 8 9 12]

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
nroi    = length(rois);

% parameters
Ntrials         = 18;
BOLDdelay       = 3; % [TR], 6 sec of assumed BOLD delay

%%%%%%%%%%%%%%%%%%%%%%%
% process starts here %
%%%%%%%%%%%%%%%%%%%%%%%
for ROI_index=1%:length(rois)% repeat this process for #ROI times
    for sub=1:Nsbj
        clear dummy voxels
        roiName = rois{ROI_index};
        Sb      = sbjs{sub};
        postfix = ['_lin_zscore_wbaseln_bd' num2str(BOLDdelay)];
        
        
        disp(['#### ' Sb ' ####']);
        disp(['++++ ' roiName ' ++++']);
        
        % Calculate number of runs: data stored in table T
        tmp         = load([condDir Sb '/nTRs.mat']);
        runID       = cellfun(@(x) contains(x,condpfx),tmp.scansdir);
        Ntrs        = unique(tmp.TRs(runID));
        Nruns       = length(tmp.TRs(runID));
        filedir     = tmp.scansdir(runID);
        cond        = struct; % to save experimental conditions
        cond.ntrial = Ntrials;
        
        for run = 1:Nruns
            clear conditionOrder* USOrder* likeOrder*
            
            %%%%
            % load data
            tmp = dir([condDir Sb '/' filepfx num2str(run) '.mat']);
            load([tmp.folder '/' tmp.name]);
            %%%%
            
            % extract run info (onsets)
            if sub==1; tseq = conditionOrder(trialOrder((run-1)*Ntrials + 1:run*Ntrials));
            else; tseq = conditionOrder_trials((run-1)*Ntrials + 1:run*Ntrials); end
            cond.TO(run,:) = round(timing_rec.trial_onset/2);
            cond.MO(run,:) = round(timing_rec.moviestart_onset/2);
            cond.S1(run,:) = round(sort([timing_rec.soundelement1_onset(tseq~=2) timing_rec.soundelement2_onset(tseq==2)])/2);
            cond.S2(run,:) = round(sort([timing_rec.soundelement1_onset(tseq==2) timing_rec.soundelement2_onset(tseq==1) timing_rec.soundelement3_onset(tseq==3)])/2);
            cond.S3(run,:) = round(sort([timing_rec.soundelement2_onset(tseq==3) timing_rec.soundelement3_onset(tseq~=3)])/2);
            cond.TS(run,:) = round(timing_rec.truckscene_onset/2);
            cond.IT(run,:) = round(timing_rec.truckscene_offset/2);
            cond.ITlen(run,:) = [timing_rec.trial_onset(2:end) - timing_rec.truckscene_offset(1:end-1) 4.2];
            
            cond.Sseq(run,:) = tseq;
            cond.Likt(run,:) = likeOrder_trials((run-1)*Ntrials + 1:run*Ntrials);
            if sub==1; cond.UStr(run,:) = USOrder(trialOrder((run-1)*Ntrials + 1:run*Ntrials));
            else; cond.UStr(run,:) = USOrder_trials((run-1)*Ntrials + 1:run*Ntrials); end
            
            cond.TRtr(run,:) = round(timing_rec.trial_onset + timing_rec.trial_offset)/2;
        end
        
        
        % load slice time course for ROI
        fprintf('loading slice time course from ROI ...');
        load([patternDir Sb '/' roiName '_fulldata.mat'],'rdata');
        rdata = rdata(idx);
        Nvoxels = size(rdata{1},1);
        if Nvoxels > 1
                       
            
            %%% preprocessing %%%
            % remove trend and normalize data
            rmrdata        = cell(size(rdata));
            nrmrdata       = cell(size(rdata));
            anrmrdata      = cell(size(rdata,1),5); % truck scene, sounds (3), whole sound sequence
            baselineMeanAll = []; baselineSTDAll = [];
            
            disp('averaging for each stimulus duration ...');
            for RUN=1:Nruns
                for STIM=1:Ntrials
                    if STIM==1
                        % detrending
                        disp('removing linear trend ...');
                        rmrdata{RUN} = remove_trend(rdata{RUN},1,0,0);
                        % normalisation
                        disp('voxel intensity normalisation')
                        baselineData = rmrdata{RUN};
                        baselineMean = mean(baselineData,2);
                        baselineSTD = std(baselineData,[],2);
                        for VOXEL=1:Nvoxels
                            nrmrdata{RUN}(VOXEL,:) = (rmrdata{RUN}(VOXEL,:) - baselineMean(VOXEL))./baselineSTD(VOXEL);
                        end
                    end
                    %%% selection of data points for decoding analyses %%%
                    anrmrdata{RUN,1}(:,STIM) = mean(nrmrdata{RUN}(:,cond.TS(RUN,STIM)+BOLDdelay:cond.TS(RUN,STIM)+2+BOLDdelay),2); % truck scene data
                    anrmrdata{RUN,2}(:,STIM) = mean(nrmrdata{RUN}(:,cond.S1(RUN,STIM)+BOLDdelay:cond.S1(RUN,STIM)+1+BOLDdelay),2); % sound 1 data
                    anrmrdata{RUN,3}(:,STIM) = mean(nrmrdata{RUN}(:,cond.S2(RUN,STIM)+BOLDdelay:cond.S2(RUN,STIM)+1+BOLDdelay),2); % sound 2 data
                    anrmrdata{RUN,4}(:,STIM) = mean(nrmrdata{RUN}(:,cond.S3(RUN,STIM)+BOLDdelay:cond.S3(RUN,STIM)+1+BOLDdelay),2); % sound 3 data
                    anrmrdata{RUN,5}(:,STIM) = mean(nrmrdata{RUN}(:,cond.MO(RUN,STIM)+1+BOLDdelay:cond.MO(RUN,STIM)+6+BOLDdelay),2); % sound full sequence data
                    anrmrdata{RUN,6}(:,STIM) = mean(nrmrdata{RUN}(:,cond.TS(RUN,STIM)+BOLDdelay:cond.TS(RUN,STIM)+3+BOLDdelay),2); % truck scene + ITI data
                end
            end
            
            
            % mean intensity pattern in ROI. this is used during the online
            % part of the experiment to calculate the spatial correlation
            % between current alignment and previous decoder construction
            % alignment (current intensity pattern and mean pattern in decoder
            % construction)
            meanpattern = [rmrdata{:}];
            meanpattern = mean(meanpattern,2);
            %Check if directories to save data exist, if not, create
            if ~exist([patternDir Sb], 'dir')
                mkdir([patternDir Sb]);
            end
            
            % remove previous variable "data"
            clear data
            disp('sorting data based on task conditions and responses...');
            sanrmrdata = sort_data( anrmrdata, cond );
            
            % info, data
            disp('saving preprocessed fMRI data ...');
            data = sanrmrdata;
            
            % Save data and mean pattern
            save([patternDir Sb '/rin_' roiName '_' condpfx(1:3) '_fMRI' postfix],'data');
            save([patternDir Sb '/rin_' roiName '_' condpfx(1:3) '_fMRI' postfix '_meanpattern'],'meanpattern');
            
        end
    end
end
