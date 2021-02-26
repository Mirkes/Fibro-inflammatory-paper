%% Load data
load('Total.mat')
% Standardise data
TotalData = zscore(TotalData);
map = rect2DMap(30, 30);
init(map, TotalData, 'pci');
%% Calculate loadings of three principal components
[~, D, V] = svds(TotalData, 3);
D = (sort(diag(D),'descend') .^ 2) / (size(TotalData, 2) - 1) / size(TotalData, 1);
%%
drawMap(map, TotalData, 'lineWidth', 0, 'labels', TotalIDs);
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);
%%
saveFigures('figures/TotalTotal.png');
%%
colours = ['b', 'g', 'r', 'k'];
lab = TotalGroups;
lab2 = lab;
lab2(lab2>1) = 1;
%%
drawMap(map, TotalData, 'classes', lab2, 'markColour', colours,  'lineWidth', 0, 'labels', TotalIDs);
legend('Healthy Volunteers (controls)', 'T2D (cases)', 'Location', 'north');
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);
%%
saveFigures('figures/Totalbinary.png');
%%
drawMap(map, TotalData, 'classes', lab, 'markColour', colours,  'lineWidth', 0, 'labels', TotalIDs);
legend('Healthy Volunteers (controls)', 'T2D standard care', 'T2D exercise', 'T2D MRP', 'Location', 'north');
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);
%%
saveFigures('figures/Totalfour.png');
