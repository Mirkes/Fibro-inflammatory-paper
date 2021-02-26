%% Load data
load('Bio.mat')
% Standardise data
BioData = zscore(BioData);
% Create map
map = rect2DMap(30, 30);
init(map, BioData, 'pci');
%% Calculate loadings of three principal components
[~, D, V] = svds(BioData, 3);
D = (sort(diag(D),'descend') .^ 2) / (size(BioData, 2) - 1) / size(BioData, 1);
%%
drawMap(map, BioData, 'lineWidth', 0, 'labels', BioIDs);
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);

%%
saveFigures('figures/BioTotal.png');
%%
colours = ['b', 'g', 'r', 'k'];
lab = BioGroups;
lab2 = lab;
lab2(lab2>1) = 1;
%%
drawMap(map, BioData, 'classes', lab2, 'markColour', colours,  'lineWidth', 0, 'labels', BioIDs);
legend('Healthy Volunteers (controls)', 'T2D (cases)', 'Location', 'north');
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);
%%
saveFigures('figures/Biobinary.png');
%%
drawMap(map, BioData, 'classes', lab, 'markColour', colours,  'lineWidth', 0, 'labels', BioIDs);
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);
%%
saveFigures('figures/Biofour.png');

