%% Load data Select one of two possible sources
load('Changes.mat');
%load('ChangesFiltered.mat');

%% Complete data in each matrix separately
[dataO, uncO] = kNNImpute(Original, 5);
[dataA, uncA] = kNNImpute(After, 5);

%% Form four different datasets
% Dataset to train maps and form normalisation
before = dataO(:, 2:end);
afters = dataA(:, 2:end);
ver = 2;

if ver == 1
    % The first version - before data used for normalisation and PC calculation only
    % Form normalisation coefficients
    means = mean(before);
    stds = std(before);
end
if ver == 2
    % The second version: all data used fro normalisation and PC calculation
    % Form normalisation coefficients
    means = mean([before; afters]);
    stds = std([before; afters]);
end

% Normalise all data
before = (before - means) ./ stds;
afters = (afters - means) ./ stds;

if ver == 1
    % The first version - before data used for normalisation and PC calculation only
    toTrain = before;
end
if ver == 2
    % The second version: all data used fro normalisation and PC calculation
    toTrain = [before; afters];
end

% Form indices of three groups
% Control
ind = Original(:, 1) == 2;
control = [before(ind, :); afters(ind, :)];
controlIDs = cellstr(['B ' + Ids(ind); 'A ' + Ids(ind)]);
ind = sum(ind);
contrClass = [ones(ind, 1); 2 * ones(ind, 1)];
% Exercise
ind = Original(:, 1) == 3;
exercise = [before(ind, :); afters(ind, :)];
exerIDs = cellstr(['B ' + Ids(ind); 'A ' + Ids(ind)]);
ind = sum(ind);
exerClass = [ones(ind, 1); 2 * ones(ind, 1)];
% MRP
ind = Original(:, 1) == 4;
mrp = [before(ind, :); afters(ind, :)];
mrpIDs = cellstr(['B ' + Ids(ind); 'A ' + Ids(ind)]);
ind = sum(ind);
mrpClass = [ones(ind, 1); 2 * ones(ind, 1)];
% Prepare colours for figures
colours = ['b', 'r', 'g', 'k'];

%% Calculate loadings of three principal components
[~, D, V] = svds(toTrain, 3);
D = (sort(diag(D),'descend') .^ 2) / (size(toTrain, 2) - 1) / size(toTrain, 1);

%% Form the map
map = rect2DMap(30, 30);
init(map, toTrain, 'pci');

%% Draw the first figure
drawMap(map, control, 'classes', contrClass, 'markColour', colours,  'lineWidth', 0, 'labels', controlIDs);
view([11.55, 2.62]);
legend('After', 'Before', 'Location','northeast');
title('T2D Standard care');
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);

%% Save figures
saveFigures('Figures/Control.png');
savefig(gcf, 'Figures/Control.fig', 'compact');

%% Draw the second figure
drawMap(map, exercise, 'classes', exerClass, 'markColour', colours,  'lineWidth', 0, 'labels', exerIDs);
view([11.55, 2.62]);
legend('After', 'Before', 'Location','northeast');
title('T2D Exercise');
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);

%% Save figures
saveFigures('figures/Exercise.png');

%% Draw the Third figure
drawMap(map, mrp, 'classes', mrpClass, 'markColour', colours,  'lineWidth', 0, 'labels', mrpIDs);
view([11.55, 2.62]);
legend('After', 'Before', 'Location','northeast');
title('T2D MRP');
xlabel(['PC1 (', sprintf('%6.4f',D(1)),')']);
ylabel(['PC2 (', sprintf('%6.4f',D(2)),')']);
zlabel(['PC3 (', sprintf('%6.4f',D(3)),')']);

%% Save figures
saveFigures('figures/MRP.png');
