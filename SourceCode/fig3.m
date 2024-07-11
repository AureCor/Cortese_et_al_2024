clearvars;

% subj info
nsb     = 41;

% rois
rois    = {'HPC', 'dlPFC'};
nr      = length(rois);

% session
session = {'Acquisition' 'Short-term test' 'Long-term test'};
ns      = length(session);

data    = csvread('../SourceData/fig3b.csv',1,2);

% data labels: "hpc_acq_Sp" "dlpfc_acq_Sp" "hpc_ext_Sp" "dlpfc_ext_Sp" "hpc_tes_Sp" "dlpfc_tes_Sp"...
%     "hpc_acq_Up" "dlpfc_acq_Up" "hpc_ext_Up" "dlpfc_ext_Up" "hpc_tes_Up" "dlpfc_tes_Up"


%% plot fig3b

% open figure

fh = figure; box off; hold on;
set(fh, 'Position', [680 550 1100 450], 'Color', 'w');
nsubplot = 2;
rix = [-.2 .2];

clr = [0 90 181;
    220 50 32]./255;

% loop through the plots/data
for i = 1:nsubplot
    subplot(1,2,i); hold on;
    
    plot([.3 3.7],[.5 .5],'k-.');
    
    for ss = 1:ns
        for R = 1:length(rois)
            
            % format data for plotting and stat test
            y = data(:,R + 2*(ss-1) + 6*(i-1));
            
            b(R) = bar(ss+rix(R), nanmean(y), .25, 'FaceColor', clr(R,:), 'EdgeColor', 'none'); alpha(b(R), .35)
            scatter(ss+rix(R) + randn(length(y),1)*.05, y, 20, clr(R,:))
            errorbar(ss+rix(R), nanmean(y), nanstd(y)/sqrt(sum(~isnan(y))), 'k', 'CapSize', 0, 'LineStyle', 'none')
            
            % test against chance, each ROI and each session
            [p(R),~,z(R)] = signrank(y, 0.5, 'tail', 'right');
        end
        % test between ROIs, each session
        y = [data(:,1 + 2*(ss-1) + 6*(i-1)) data(:,2 + 2*(ss-1) + 6*(i-1))];
        [p(3),~,z(3)] = signrank(y(:,1),y(:,2));
        [~, ~, ~, adj_p] = fdr_bh(p);
        disp([rois{1} ' [' session{ss} '] z=' num2str(z(1).zval,3) ', p=' num2str(p(1),3) ', multcmp BH p=' num2str(adj_p(1),2)])
        disp([rois{2} ' [' session{ss} '] z=' num2str(z(2).zval,3) ', p=' num2str(p(2),3) ', multcmp BH p=' num2str(adj_p(2),2)])
        disp(['Between ROIs,' session{ss} ']: z=' num2str(z(3).zval,3) ', p=' num2str(p(3),3) ', multcmp BH p=' num2str(adj_p(3),2)])
    end
    
    
    
    if i==1; legend(b, rois, 'box', 'off', 'location', 'northwest'); end
    % set axes, limits, ticks, labels, etc
    ax = gca;
    set(ax, 'ylim', [.2 .9], 'ytick', round(.25:.125:.9,2), 'xtick', 1:3, 'xticklabels', {'Acquisition' 'Short-term\newlinetest' 'Long-term\newlinetest'});
    ax.YLabel.String = 'decoding accuracy';
    if i==1; text(.7,.93,'Sound sequence period','FontSize',14); else; text(0.8,.93,'US anticipation epoch','FontSize',14); end
end

%% Linear mixed effect models
% -------------------------------------------------------------------------

%% [S-cs+ vs E-cs+, sound sequenced & US period]

y = data(:)-.5;
sub = (repmat((1:nsb)',nr*ns*2,1));
roi = categorical(repmat(reshape(repmat((1:nr),nsb,1),nsb*nr,1),ns*2,1));
ses = categorical(repmat(reshape(repmat((1:ns),nsb*nr,1),nsb*nr*ns,1),2,1));
epo = categorical([ones(nsb*nr*ns,1); 2*ones(nsb*nr*ns,1)]);

t = table(y, sub, roi, ses, epo);

f = 'y ~ roi*ses*epo + (roi*ses*epo | sub)';
glme1 = fitglme(t,f);

f = 'y ~ roi*ses*epo - roi:ses:epo + (roi*ses*epo - roi:ses:epo | sub)';
glme2 = fitglme(t,f);

cg = compare(glme2, glme1);
if cg.pValue<.05; disp(glme1); model = glme1; else; disp(glme2); model = glme2; end


%% [S-cs+ vs E-cs+, sound sequence]
f = 'y ~ roi*ses + (roi*ses | sub)';
glme1 = fitglme(t(t.epo=='1',:),f);

f = 'y ~ roi+ses + (roi+ses | sub)';
glme2 = fitglme(t(t.epo=='1',:),f);

cg = compare(glme2, glme1);
if cg.pValue<.05; disp(glme1); model = glme1; else; disp(glme2); model = glme2; end


%% [S-cs+ vs E-cs+, US period]
f = 'y ~ roi*ses + (roi*ses | sub)';
glme1 = fitglme(t(t.epo=='2',:),f);

f = 'y ~ roi+ses + (roi+ses | sub)';
glme2 = fitglme(t(t.epo=='2',:),f);

cg = compare(glme2, glme1);
if cg.pValue<.05; disp(glme1); model = glme1; else; disp(glme2); model = glme2; end
