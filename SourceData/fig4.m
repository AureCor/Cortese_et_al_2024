% script to plot the results of the sparse linear regression to evaluate
% functional coupling at the multixovel pattern level between two regions

clearvars;

% subj info
nsb      = 41;

sessions = {'Acquisition' 'Immediate\newlinetest' 'Long-term\newlinetest'};

sROI     = {'HPC' 'BA46_9'}; % source ROI (likelihood)
tROI     = {'amyg-vmPFC'}; % target ROI (timecourse)

sROI_tit = {'HPC', 'DLPFC'}; % source ROI (likelihood) - title

% load data
tmp      = readtable('../SourceData/fig4b.csv');
data     = table2array(tmp(:,3:end));
% data labels: "hpc_acq" "dlpfc_acq" "hpc_ext" "dlpfc_ext" "hpc_tes"
% "dlpfc_tes"



%% Plot fig4b
disp('--------------------')

% plot params
idx_r = [1 3 5;
         2 4 6];
idx_x = [1 2 3];
clr = [0 90 181;
       220 50 32]./255;

% open figure space
fh = figure('Position',[200 200 800 600],'color','w'); hold on
plot([.5 3.5], [0 0], 'k-.')
    
for ss = 1:length(sessions)
    for R=1:length(sROI_tit)

        % format data for plotting and stat test
        y = data(:,R + 2*(ss-1));

        scatter(ss+.2*(R-1), mean(y), 150, clr(R,:), 's', 'filled'); 
        s = scatter(ss+.2*(R-1)+randn(length(y),1)*.05, y, 40, clr(R,:), 's', 'filled'); alpha(s, .2); 
        % test against chance, each ROI and each session
        [p(R),~,z(R)] = signrank(y, 0);

    end
    % test between ROIs, each session
    y = [data(:,1 + 2*(ss-1)) data(:,2 + 2*(ss-1))];
    [p(3),~,z(3)] = signrank(y(:,1),y(:,2));
    [~, ~, ~, adj_p] = fdr_bh(p);
    
    disp([sROI_tit{1} ' [' sessions{ss} '] z=' num2str(z(1).zval,3) ', p=' num2str(p(1),3) ', multcmp BH p=' num2str(adj_p(1),2) ' , r=' num2str(abs(z(1).zval)/sqrt(nsb),2)])
    disp([sROI_tit{2} ' [' sessions{ss} '] z=' num2str(z(2).zval,3) ', p=' num2str(p(2),3) ', multcmp BH p=' num2str(adj_p(2),2) ' , r=' num2str(abs(z(2).zval)/sqrt(nsb),2)])
    disp(['Between ROIs,' sessions{ss} ']: z=' num2str(z(3).zval,3) ', p=' num2str(p(3),3) ', multcmp BH p=' num2str(adj_p(3),2) ' , r=' num2str(abs(z(3).zval)/sqrt(nsb),2)])
    clear z p adj_p

end

for R=1:length(sROI_tit)
    y = data(:,idx_r(R,:));
    e(R) = errorbar(idx_x+.2*(R-1), mean(y), std(y)/sqrt(nsb),...
        'Color',clr(R,:), 'Linestyle', '-','LineWidth',2,'Capsize',0);
    
    [p(1),~,z(1)] = signrank(y(:,1), y(:,2));
    [p(2),~,z(2)] = signrank(y(:,1), y(:,3));
    [p(3),~,z(3)] = signrank(y(:,2), y(:,3));
    [~, ~, ~, adj_p] = fdr_bh(p);
    disp([sROI_tit{R} ' [acquisition vs immediate]: z=' num2str(z(1).zval,3) ', p=' num2str(p(1),3) ', multcmp BH p=' num2str(adj_p(1),2) ' , r=' num2str(abs(z(1).zval)/sqrt(nsb),2)])
    disp([sROI_tit{R} ' [acquisition vs long-term]: z=' num2str(z(2).zval,3) ', p=' num2str(p(2),3) ', multcmp BH p=' num2str(adj_p(2),2) ' , r=' num2str(abs(z(2).zval)/sqrt(nsb),2)])
    disp([sROI_tit{R} ' [immediate vs long-term]: z=' num2str(z(3).zval,3) ', p=' num2str(p(3),3) ', multcmp BH p=' num2str(adj_p(3),2) ' , r=' num2str(abs(z(3).zval)/sqrt(nsb),2)])
    clear z p adj_p
end

xlim([0 length(sessions)+1]);

set(gca, 'XTickLabel',sessions(idx_x),'XTick',idx_x+.1, 'ytick', -1:.25:1.5, 'FontSize', 12);
ylabel('Fisher-transformed \newlinecorrelation coefficient'); 
legend(e, sROI_tit, 'Box', 'off', 'Location', 'best');



%% Linear mixed effect model
% -------------------------------------------------------------------------

disp('--------------------')

idx  = [1 2 3]; ns = 3; nr = 2;

y    = data(:);
sub  = repmat((1:nsb)',ns*nr,1);
time = categorical([ones(nsb*nr,1); 2*ones(nsb*nr,1); 3*ones(nsb*nr,1)]);
roi  = categorical(repmat([ones(nsb,1); 2*ones(nsb,1)],ns,1));

%load decoding accuracy data
tmp  = readtable('../SourceData/fig4b_dac.csv');
dac  = table2array(tmp);
dac  = dac(:);

t = table(y,sub,time,roi,dac);

f = 'y ~ roi*time + (roi*time|sub) + (1|dac)';
glme1 = fitglme(t,f);
disp(glme1)


%% Anxiety vs diff in ealk HPC-DLPFC long-term test 

% load data
data = readtable('../SourceData/fig4c.csv');
x1   = data.vmpfc;
x2   = data.amyg;
y    = data.traitAnxiety;

disp('----------------------------------')
fh = figure('Position',[200 200 600 500],'color','w'); hold on;

plot(x1, y, 'k^')
plot(x2, y, 'ro')
lsline; 
ax = gca; ax.XLabel.String = '\Delta corr [DLPFC-HPC] at long-term test';
ax.YLabel.String = 'Trait anxiety';
ax.XLim = [-1 2];

legend({'vmPFC' 'Amyg'}, 'box', 'off')

disp('vmPFC')
[r1,p] = corr(x1,y,'type','Pearson');
disp(['Pearson r = ' num2str(r1,3) ', p = ' num2str(p,3)])
[a,b] = robustfit(x1,y);
disp(['Robust regression slope = ' num2str(a(2),3) ', t(' num2str(b.dfe) ') = ' num2str(b.t(2),3) ', p = ' num2str(b.p(2),3)])

disp('Amyg')
[r2,p] = corr(x2,y,'type','Pearson');
disp(['Pearson r = ' num2str(r2,3) ', p = ' num2str(p,3)])
[a,b] = robustfit(x2,y);
disp(['Robust regression slope = ' num2str(a(2),3) ', t(' num2str(b.dfe) ') = ' num2str(b.t(2),3) ', p = ' num2str(b.p(2),3)])

[p, z, za, zb] = corr_rtest(r1, r2, 41, 41);
disp('comparison')
disp(['z = ' num2str(z,3) ', p = ' num2str(p,3)])
