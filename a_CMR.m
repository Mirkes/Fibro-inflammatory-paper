%% Load data
load('CMR.mat')
% Standardise data
CMRData = zscore(CMRData);
map = rect2DMap(30, 30);
init(map, CMRData, 'pci');
%% Calculate loadings of three principal components
[~, D, V] = svds(CMRData, 3);
D = (sort(diag(D),'descend') .^ 2) / (size(CMRData, 2) - 1) / size(CMRData, 1);
%%
drawMap(map, CMRData, 'lineWidth', 0, 'labels', CMRIDs);
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);
%%
saveFigures('figures/CMRTotal.png');
%%
colours = ['b', 'g', 'r', 'k'];
lab = CMRGroups;
lab2 = lab;
lab2(lab2>1) = 1;
%%
drawMap(map, CMRData, 'classes', lab2, 'markColour', colours,  'lineWidth', 0, 'labels', CMRIDs);
legend('Healthy Volunteers (controls)', 'T2D (cases)', 'Location', 'north');
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);
%%
saveFigures('figures/CMRbinary.png');
%%
drawMap(map, CMRData, 'classes', lab, 'markColour', colours,  'lineWidth', 0, 'labels', CMRIDs);
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);
%%
saveFigures('figures/CMRfour.png');
