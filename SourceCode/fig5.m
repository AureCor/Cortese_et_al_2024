% script to plot the results of the sparse linear regression to evaluate
% functional coupling at the multixovel pattern level between two regions

clearvars;

% subj info
nsb    = 41;

sessions = {'Acquisition' 'Immediate\newlinetest' 'Long-term\newlinetest'};

sROI = {'HPC' 'BA46_9'}; % source ROI (likelihood)
tROI = {'vmPFC' 'amyg'}; % target ROI (timecourse)
trc  = 1;

sROI_tit = {'HPC-vmPFC', 'DLPFC-vmPFC'}; % source ROI (likelihood) - title

% load data
tmp      = readtable('../SourceData/fig5.csv');
data     = table2array(tmp(:,3:end));
% data labels: "cspe_hpc_acq" "cspe_dlpfc_acq" "cspe_hpc_ext"
% "cspe_dlpfc_ext" "cspe_hpc_tes" "cspe_dlpfc_tes" "csps_hpc_acq" 
% "csps_dlpfc_acq" "csps_hpc_ext" "csps_dlpfc_ext" "csps_hpc_tes" "csps_dlpfc_tes"

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------')
% plot params
clr = [0 90 181;
    220 50 32]./255;
idx_r = [1 3 5;
         2 4 6];
idx_x = [1 2 3];

for jj = 1:2 % there are two control analyses, 1st is CS+seq vs CS-, 2nd is CS+ele vs CS-
    
    % open figure space
    fh = figure('Position',[680 550 800 600],'color','w'); hold on;
    plot([.5 3.5], [0 0], 'k-.')
    
    for ss = 1:length(sessions)
        
        for R=1:length(sROI_tit)

            % format data for plotting and stat test
            y = data(:,R + 2*(ss-1) + 6*(jj-1));

            scatter(ss+.2*(R-1), mean(y), 150, clr(R,:), 'o', 'filled');
            s = scatter(ss+.2*(R-1)+randn(length(y),1)*.05, y, 40, clr(R,:), 'o', 'filled'); alpha(s, .2);
            % test against chance, each ROI and each session
            [p(R),~,z(R)] = signrank(y, 0);

        end

        y = [data(:,1 + 2*(ss-1) + 6*(jj-1)) data(:,2 + 2*(ss-1) + 6*(jj-1))];
        [p(3),~,z(3)] = signrank(y(:,1),y(:,2));
        [~, ~, ~, adj_p] = fdr_bh(p);
        
        disp([sROI_tit{1} ' [' sessions{ss} '] z=' num2str(z(1).zval,3) ', p=' num2str(p(1),3) ', multcmp BH p=' num2str(adj_p(1),2)])
        disp([sROI_tit{2} ' [' sessions{ss} '] z=' num2str(z(2).zval,3) ', p=' num2str(p(2),3) ', multcmp BH p=' num2str(adj_p(2),2)])
        disp(['Between ROIs,' sessions{ss} ']: z=' num2str(z(3).zval,3) ', p=' num2str(p(3),3) ', multcmp BH p=' num2str(adj_p(3),2)])
        clear z p adj_p

    end
    
    for R=1:length(sROI_tit)

        y = data(:,idx_r(R,:) + 6*(jj-1));

        e(R) = errorbar(idx_x+.2*(R-1), mean(y), std(y)/sqrt(nsb),...
            'Color',clr(R,:), 'Linestyle', '-','LineWidth',2,'Capsize',0);
        
        [p(1),~,z(1)] = signrank(y(:,1), y(:,2));
        [p(2),~,z(2)] = signrank(y(:,1), y(:,3));
        [p(3),~,z(3)] = signrank(y(:,2), y(:,3));
        [~, ~, ~, adj_p] = fdr_bh(p);
        disp([sROI_tit{R} ' [acquisition vs immediate]: z=' num2str(z(1).zval,3) ', p=' num2str(p(1),3) ', multcmp BH p=' num2str(adj_p(1),2)])
        disp([sROI_tit{R} ' [acquisition vs long-term]: z=' num2str(z(2).zval,3) ', p=' num2str(p(2),3) ', multcmp BH p=' num2str(adj_p(2),2)])
        disp([sROI_tit{R} ' [immediate vs long-term]: z=' num2str(z(3).zval,3) ', p=' num2str(p(3),3) ', multcmp BH p=' num2str(adj_p(3),2)])
        clear z p adj_p

    end
    
    xlim([0 length(sessions)+1]);
    
    set(gca, 'XTickLabel',sessions(idx_x),'XTick',idx_x+.1, 'ytick', -1:.25:1.5, 'ylim', [-1 1.5], 'FontSize', 12);
    ylabel('Fisher-transformed \newlinecorrelation coefficient');
    legend(e, sROI_tit, 'Box', 'off', 'Location', 'best');
    
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------')
idxn = {'b' 'd'};
idx = [1 2 3]; ns = 3; nr = 2;

%load decoding accuracy data
tmp  = readtable('../SourceData/fig5_dac.csv');
dac  = table2array(tmp);
dac  = dac(:);

for jj = 1:2 % 1st is CS+seq vs CS-, 2nd is CS+ele vs CS-
    
    y = data(:,[1 2 3 4 5 6] + 6*(jj-1));
    y   = y(:);
    sub = repmat((1:nsb)',ns*nr,1);
    time = categorical([ones(nsb*nr,1); 2*ones(nsb*nr,1); 3*ones(nsb*nr,1)]);
    roi = categorical(repmat([ones(nsb,1); 2*ones(nsb,1)],ns,1));
    
    t = table(y,sub,time,roi,dac);
    
    f = 'y ~ roi*time + (roi*time|sub) + (1|dac)';
    glme = fitglme(t,f);
    disp(glme)

end