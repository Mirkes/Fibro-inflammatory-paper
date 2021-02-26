% Search of the best separation of two groups
%
% Data preparation
%
%% Load data Select one of two possible sources
%load('Changes.mat');
load('ChangesFiltered.mat');

%% Complete data in each matrix separately
[dataO, ~] = kNNImpute(Original, 5);
[dataA, ~] = kNNImpute(After, 5);

%% Form four different datasets
% Dataset to use
before = dataO(:, 2:end);
afters = dataA(:, 2:end);
Columns = Columns(2:end);

%% Form set of three pairs of datasets
dbNames = {'Standard care', 'Excersise', 'MRP'};
dbs = cell(3, 3);
for group = 1:3
    ind = After(:, 1) == group + 1;
    dbs(group, 1) = dbNames(group);
    dbs(group, 2) = {before(ind, :)};
    dbs(group, 3) = {afters(ind, :)};
end
clear After dataA dataO dbNames ids ind k Original

%% Search of the best one attribute classifier for all three groups
N = length(Columns);
for group = 1:3
    % extract data from cells
    name = dbs{group, 1};
    class1 = dbs{group, 2};
    class2 = dbs{group, 3};
    bestAtr = -1;
    bestErr = inf;
    bestThr = -1;
    for kk = 1:N
        [thr, err] = oneDClass(class1(:, kk), class2(:, kk));
        if err < bestErr
            bestErr = err;
            bestThr = thr;
            bestAttr = kk;
        end
    end
    fprintf(['For group\t"%s"\tthe best attribute is\t%d\t"%s"\twith',...
        ' threshold\t%f\tand error\t%f\n'],...
        name, bestAttr, Columns(bestAttr), bestThr, bestErr);
    [thr, err] = oneDClass(class1(:, bestAttr), class2(:, bestAttr), [Columns(bestAttr), 'Before', 'After']);
    saveFigures(['figures/1D_',name,'.png']);
    close();
end

%% Search the best 2D Fisher discriminant separation
N = length(Columns);
mnemo = {'STD', 'EXC', 'MRP'};
for group = 1:3
    % extract data from cells
    name = dbs{group, 1};
    class1 = dbs{group, 2};
    class2 = dbs{group, 3};
    bestAtr = [];
    bestErr = inf;
    bestThr = -1;
    bestDir = [];
    for k = 1:N - 1
        for kk = k + 1:N
            [thr, err, dir] = fisher(class1(:, [k, kk]), class2(:, [k, kk]));
            if err < bestErr
                bestErr = err;
                bestThr = thr;
                bestAttr = [k, kk];
                bestDir = dir;
            end
        end
    end
    % Print result to console
    fprintf(['For group\t"%s"\tthe best attributes are\t%d\t"%s"\tand\t',...
        '%d\t"%s"\twith threshold\t%f\tand error\t%f\n'],...
        name, bestAttr(1), Columns(bestAttr(1)), bestAttr(2),...
        Columns(bestAttr(2)), bestThr, bestErr);
    % Draw histograms for the best couple of attributes
    [thr, err] = fisher(class1(:, bestAttr), class2(:, bestAttr),...
        [convertCharsToStrings([char(Columns(bestAttr(1))), ' and ',...
        char(Columns(bestAttr(2)))]), 'Before', 'After']);
    saveFigures(['figures/2D_',name,'.png']);
    close();
    % Draw distribution of points for the bect couple of attributes
    plot(class1(:, bestAttr(1)), class1(:, bestAttr(2)), 'sb', 'MarkerFaceColor', 'b');
    hold on;
    plot(class2(:, bestAttr(1)), class2(:, bestAttr(2)), 'sr', 'MarkerFaceColor', 'r');
    legend('Before', 'After');
    xlabel(Columns(bestAttr(1)));
    ylabel(Columns(bestAttr(2)));
    title(name);
    % Draw separation line
    % Calculate centre
    cent = bestDir' * bestThr;
    % calculate direction of separability line
    dir = [-bestDir(2), bestDir(1)];
    % Get limits of drawing area
    xl = xlim();
    yl = ylim();
    % Calculate marginal values of parameter in equation p = centr + t*dir
    t = sort([([xl(1), yl(1)] - cent) ./ dir, ([xl(2), yl(2)] - cent) ./ dir]);
    tmp = [cent + t(2) * dir; cent + t(3) * dir];
    plot(tmp(:, 1), tmp(:, 2), '-k');
    xlim(xl);
    ylim(yl);
    saveFigures(['figures/', mnemo{group}, '_plain.png']);
    close();
end

%% The full dimensional Fisher’s discriminant
for group = 1:3
    % extract data from cells
    name = dbs{group, 1};
    class1 = dbs{group, 2};
    class2 = dbs{group, 3};
    [thr, err, dir] = fisher(class1, class2, {name, 'Before', 'After'});
    fprintf(['For group\t"%s"\tthe best threshold\t%f\tand error\t%f\n'],...
        name, thr, err);
end

%% Dimension analysis
% Join all databases
db = [dbs{1, 2}; dbs{1, 3}; dbs{2, 2}; dbs{2, 3}; dbs{3, 2}; dbs{3, 3}];
means = mean(db);
stds = std(db);
% standardize each database
db = (db - means) ./ stds;
for group = 1:3
    for k = 2:3
        dbs{group, k} = (dbs{group, k} - means) ./ stds;
    end
end

%% calculation of eigenvalues
s = sort(eig(cov(db)), 'descend');
figure;
plot(s);
hold on;
plot([0, size(db, 2)],[1,1],'-k');
[~, g] = brokenStick(s);
plot(g * size(db, 2), '-r');
xlabel('The Number of principal components');
ylabel('Fraction of variance explained');
set(gca,'fontsize', 14);
legend('Observed', 'Kaiser rule', 'Brocken stick');
saveFigures('figures/PCs.png');
[~, sepDim] = SeparabilityAnalysis(db);
