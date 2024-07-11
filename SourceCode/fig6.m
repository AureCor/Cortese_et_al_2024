%% SCR vs connectivity [day 1, day 2]
clearvars;

% load data
data = readtable('../SourceData/fig6.csv');
x1   = data.DLPFC_vmPFCD1;
x2   = data.DLPFC_vmPFCD2;
y1   = data.SCRD1;
y2   = data.scrD2;

fh = figure('Color', 'w'); hold on; box off;
s(1) = scatter(x1, y1, 50, [0.366 .821 .673], 'filled');
s(2) = scatter(x2, y2, 50, [0.066 .221 .073], 'filled');
ax = gca; ax.XLabel.String = 'information transmission [z]';
ax.YLabel.String = '\Delta SCR (E-cs+ - cs-)';
title('DLPFC-vmPFC')

disp('day 1')
[r1,p] = corr(x1, y1,'type','Pearson');
disp(['r = ' num2str(r1,3) ', p = ' num2str(p,3)])
[a,b] = robustfit(x1,y1);
disp(['Robust regression slope = ' num2str(a(2),3) ', t(' num2str(b.dfe) ') = ' num2str(b.t(2),3) ', p = ' num2str(b.p(2),3)])

x = -.85:.001:1.25;
plot(x, a(1) + a(2)*x, 'Color', [0.366 .821 .673], 'LineWidth', 1);

disp('day 2')
[r2,p] = corr(x2,y2,'type','Pearson');
disp(['r = ' num2str(r2,3) ', p = ' num2str(p,3)])
[a,b] = robustfit(x2,y2);
disp(['Robust regression slope = ' num2str(a(2),3) ', t(' num2str(b.dfe) ') = ' num2str(b.t(2),3) ', p = ' num2str(b.p(2),3)])

x = -.85:.001:1.25;
plot(x, a(1) + a(2)*x, 'Color', [0.066 .221 .073], 'LineWidth', 1);

legend(s, {'Day 1', 'Day 2'}, 'location', 'northwest', 'box', 'off')

[p, z] = corr_rtest(r1, r2, 41, 41);
disp(['z test, difference in correlation: z = ' num2str(z,3) ', p = ' num2str(p(2),3)])
