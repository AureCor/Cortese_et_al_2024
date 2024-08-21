clearvars;

session = {'Acquisition' 'Short-term\newlinetest' 'Long-term\newlinetest'};
nsb     = 42;
ns      = length(session);
nc      = 2;

% load SCR data
scr = csvread('../SourceData/fig2.csv',1,2);


%% linear mixed effect model

% scr
y = [scr(:,[1 2])-scr(:,3) scr(:,[4 5])-scr(:,6) scr(:,[7 8])-scr(:,9)];
y = y(:);

% code variables
sub = (repmat((1:nsb)',nc*ns,1));
cst = categorical(repmat(reshape(repmat((1:nc),nsb,1),nsb*nc,1),ns,1));
ses = categorical([ones(nsb*nc,1); 2*ones(nsb*nc,1); 3*ones(nsb*nc,1)]);

% create table
t = table(y, sub, ses, cst);

% run model 
f = 'y ~ ses*cst + (ses*cst | sub)';
glme1 = fitglme(t, f);
disp(glme1)


%% fig2a

% open figure 
fh = figure; box off; hold on;
set(fh, 'Position', [680 550 500 400], 'Color', 'w');

% load color info
load('colors.mat')

for ss = 1:ns
    
    % format data for plotting and stat test
    y = mean((scr(:,[1 2]+3*(ss-1)) - scr(:,3+3*(ss-1))),2);
    
    % stats
    [p(ss),h,z(ss)] = signrank(y,0,'tail','right');
    disp([session{ss} ', cs+: z=' num2str(z(ss).zval,3) ', p=' num2str(p(ss),3) ' --- N=' num2str(length(y)) ' , effect size r=' num2str(abs(z(ss).zval)/sqrt(length(y)),2)])
    
    tmp = max(max(y(:)));
    if p(ss)<.001; text(ss, tmp+.06, '***', 'HorizontalAlignment', 'center','FontSize',12);
    elseif p(ss)<.01; text(ss, tmp+.06, '**', 'HorizontalAlignment', 'center','FontSize',12);
    elseif p(ss)<.05; text(ss, tmp+.06, '*', 'HorizontalAlignment', 'center','FontSize',12);
    elseif p(ss)<.1; text(ss, tmp+.06, num2str(p(ss),2), 'HorizontalAlignment', 'center','FontSize',9)
    else; text(ss, tmp+.06, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    
    % plot
    b = bar(ss, nanmean(y), .3, 'FaceColor', [.3 .3 .3], 'EdgeColor', 'none'); alpha(b, .35)
    s = scatter(ss + randn(length(y),1)*.05, y, 40, [.3 .3 .3], 'filled');
    if p(ss)<.1; alpha(s,.5); else; alpha(s,.15); end
    errorbar(ss, nanmean(y), nanstd(y)/sqrt(sum(~isnan(y))), 'k', 'CapSize', 0, 'LineStyle', 'none', 'LineWidth', 2)
    
end

[h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p);

disp(['BH-multcomp: p=' num2str(adj_p,3)])

legend(b, {'cs+'}, 'box', 'off', 'location', 'southwest');
% set axes, limits, ticks, labels, etc
ax = gca;
set(ax, 'ylim', [-.25 .25], 'ytick', -.3:.1:.3, 'xtick', 1:3, 'xticklabels', session);
ax.YLabel.String = '\DeltaSCR [mean(cs+) - cs-]'; 



%% fig 2b

% open figure 
fh = figure; box off; hold on;
set(fh, 'Position', [680 550 600 400], 'Color', 'w');


% load color info
load('colors.mat')
clix    = [1 2];
ixs     = [-.2 .2];

for ss = 1:ns
    
    % format data for plotting and stat test
    y = (scr(:,[1 2]+3*(ss-1)) - scr(:,3+3*(ss-1)));
    
    % compute stats
    [pu(1,ss),h,zu(1,ss)] = signrank(y(:,1),0,'tail','right');
    [pu(2,ss),h,zu(2,ss)] = signrank(y(:,2),0,'tail','right');
    [pc(ss),h,zc(ss)] = signrank(y(:,1),y(:,2));
    
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([pu(1,ss) pu(2,ss) pc(ss)]);
    
    disp([session{ss} ', S-cs+: z=' num2str(zu(1,ss).zval,3) ', p=' num2str(pu(1,ss),3) ' --- N=' num2str(length(y)) ' , effect size r=' num2str(abs(zu(1,ss).zval)/sqrt(length(y)),2)])
    disp([session{ss} ', E-cs+: z=' num2str(zu(2,ss).zval,3) ', p=' num2str(pu(2,ss),3) ' , effect size r=' num2str(abs(zu(2,ss).zval)/sqrt(length(y)),2)])
    disp([session{ss} ', comp: z=' num2str(zc(ss).zval,3) ', p=' num2str(pc(ss),3) ' , effect size r=' num2str(abs(zc(ss).zval)/sqrt(length(y)),2)])
     
    disp([session{ss} ', BH-multcomp, S-cs+: z=' num2str(zu(1,ss).zval,3) ', p=' num2str(adj_p(1),3) ' --- N=' num2str(length(y))])
    disp([session{ss} ', BH-multcomp, E-cs+: z=' num2str(zu(2,ss).zval,3) ', p=' num2str(adj_p(2),3)])
    disp([session{ss} ', BH-multcomp, comp: z=' num2str(zc(ss).zval,3) ', p=' num2str(adj_p(3),3)]) 
    
    % plot stats
    tmp = max(max(y(:)));
    
    % individual
    % CS+ sequence
    plot([-.275 -.125]+ss, [tmp tmp]+.02,'k-');
    if adj_p(1)<.001; text(ss-.2, tmp+.0275, '***', 'HorizontalAlignment', 'center', 'FontSize', 13);
    elseif adj_p(1)<.01; text(ss-.2, tmp+.0275, '**', 'HorizontalAlignment', 'center', 'FontSize', 13);
    elseif adj_p(1)<.05; text(ss-.2, tmp+.0275, '*', 'HorizontalAlignment', 'center', 'FontSize', 13);
    elseif adj_p(1)<.1; text(ss-.2, tmp+.035, num2str(adj_p(1),2), 'HorizontalAlignment', 'center', 'FontSize', 9);
    else; text(ss-.2, tmp+.035, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % CS+ element
    plot([.125 .275]+ss, [tmp tmp]+.02,'k-');
    if adj_p(2)<.001; text(ss+.2, tmp+.0275, '***', 'HorizontalAlignment', 'center', 'FontSize', 13);
    elseif adj_p(2)<.01; text(ss+.2, tmp+.0275, '**', 'HorizontalAlignment', 'center', 'FontSize', 13);
    elseif adj_p(2)<.05; text(ss+.2, tmp+.0275, '*', 'HorizontalAlignment', 'center', 'FontSize', 13);
    elseif adj_p(2)<.1; text(ss+.2, tmp+.035, num2str(adj_p(2),2), 'HorizontalAlignment', 'center', 'FontSize', 9);
    else; text(ss+.2, tmp+.035, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    % comparison
    plot([-.2 .2]+ss, [tmp tmp]+.07,'k-');
    if adj_p(3)<.001; text(ss, tmp+.0775, '***', 'HorizontalAlignment', 'center', 'FontSize', 13);
    elseif adj_p(3)<.01; text(ss, tmp+.0775, '**', 'HorizontalAlignment', 'center', 'FontSize', 13);
    elseif adj_p(3)<.05; text(ss, tmp+.0775, '*', 'HorizontalAlignment', 'center', 'FontSize', 13);
    elseif adj_p(3)<.1; text(ss, tmp+.085, num2str(adj_p(3),2), 'HorizontalAlignment', 'center', 'FontSize', 9);
    else; text(ss, tmp+.085, 'ns', 'HorizontalAlignment', 'center', 'FontSize', 10);
    end
    
    % plot data
    b(1) = bar(ss+ixs(1), nanmean(y(:,1)), .3, 'FaceColor', clr(clix(1),:), 'EdgeColor', 'none'); alpha(b(1), .35)
    b(2) = bar(ss+ixs(2), nanmean(y(:,2)), .3, 'FaceColor', clr(clix(2),:), 'EdgeColor', 'none'); alpha(b(2), .35)
    s(1) = scatter(ss+ixs(1) + randn(length(y),1)*.05, y(:,1), 40, clr(clix(1),:), 'filled'); 
    if pu(1,ss)<.1; alpha(s(1),.5); else; alpha(s(1),.15); end
    s(2) = scatter(ss+ixs(2) + randn(length(y),1)*.05, y(:,2), 40, clr(clix(2),:), 'filled'); 
    if pu(2,ss)<.1; alpha(s(2),.5); else; alpha(s(2),.15); end
    errorbar(ss+ixs(1), nanmean(y(:,1)), nanstd(y(:,1))/sqrt(sum(~isnan(y(:,1)))), 'k', 'CapSize', 0, 'LineStyle', 'none', 'LineWidth', 2)
    errorbar(ss+ixs(2), nanmean(y(:,2)), nanstd(y(:,2))/sqrt(sum(~isnan(y(:,2)))), 'k', 'CapSize', 0, 'LineStyle', 'none', 'LineWidth', 2)
    
end

legend(b, {'S-cs+' 'E-cs+'}, 'box', 'off', 'location', 'southwest');
% set axes, limits, ticks, labels, etc
ax = gca;
set(ax, 'ylim', [-.25 .28], 'ytick', -.3:.1:.3, 'xtick', 1:3, 'xticklabels', session);
ax.YLabel.String = '\DeltaSCR (csX - cs-)';