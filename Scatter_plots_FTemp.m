% scatter_plots.m
% In this script, scatter plots are generated for the pair
%  - E_{π} and T_{-0.1441}% Save plots to images? Set to true if yes
% ... with their respective regression lines.

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console
pkg load statistics;

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 2; % the line thickness in the graph
fontSize = 16; % the font size in the graph
saveToFile = true; % Set to true to auto-save plots


saveToFile = false;

% Utility: Below is a function to round-off to 4 decimal places | returns string
as_4_dp_str = @(x) sprintf('%.06f', round(x*(10^6))/(10^6));

% Cell containing Entropy and Heat Capacity of 30 lower benzenoids
expData = {reshape([% Pi Electroni Energy
   8.0000  13.6832 19.3137 19.4483 24.9308 25.1875 25.1012 25.1922 25.2745 22.5055
   30.5440 30.7255 30.8805 30.8795 30.7627 30.9990 30.9362 30.9386 30.9432 30.8390
   30.9418 30.8336 28.2453 28.3361 28.2220 36.6814 31.4251 36.1557 34.5718 46.4974
   ]', 30, 1),
   "E_{π}" % Their labels
 };

% Number of edges in each dudv partition: (2,2), (2,3), then (3,3)
d_f = [
   12 12 12 14 12 16 14 16 18 12 12 14 16 16 14 20 18 16 18 16 18 16 16 16 14 12 14 12 12 8
   0  4  8  6  12 8  10 8  6  8  16 14 12 12 14 8  10 12 10 12 10 10 8  8  10 12 10 20 12 16
   0  2  4  6  6  10 8  10 12 10 8  10 12 12 10 16 14 12 14 12 14 12 16 16 14 24 20 10 24 38
 ]';

% Define nu_exp here before using it in getIndexFns Lower 30 Bhs (nu_exp: number of vertices)
nu_exp = [
    6 10 14 14 18 18 18 18 18 18 22 22 22 22 22 22 22 22 22 22 22 22 20 20 20 24 22 26 24 32
]';

% Index-computing function for T_{a}
getIndexFns = {
    @(beta) d_f(:,1) .* (2 ./ (nu_exp - 2)).^beta + ...
             d_f(:,2) .* (2 ./ (nu_exp - 2) + 3 ./ (nu_exp - 3)).^beta + ...
             d_f(:,3) .* (3 ./ (nu_exp - 3)).^beta;
}';

% Cell containing their respective labels
indexName = {"T"};

% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % one (E_{π})
numIndices = size(getIndexFns,2); % two (T_{a}^1 & T_{a}^2)
numCurves = numData*numIndices;   % two (2 curves)

for edn = 1:numData % edn = experimental data number | 1=E_{π}
  for fnn = 1:numIndices % fnn = function number | n=1:T_{a}^1, n=2:T_{a}^2
    ccFn = @(alpha) corrcoef( % Gets corrcoef between bp and index
      getIndexFns{fnn}(alpha)(!isnan(expData{1,edn})),
      expData{1,edn}(!isnan(expData{1,edn}))
    )(1,2);

    peakAlpha = mean(
      GoldenSectionSearch_Maximum(ccFn, -5, 5, 1e-15));

    % Regression line i.e., y = mx + b;
    model = polyfit(getIndexFns{fnn}(peakAlpha), expData{1,edn},1);
    m = model(1); b = model(2);         % For the regression line
    x = [getIndexFns{fnn}(peakAlpha)(1) max(getIndexFns{fnn}(peakAlpha))];
    y = m*x + b;

    % Scatter plot
    this_figure = figure(fnn); hold on;
    regLine = plot(x, y, '-', 'LineWidth',lineWidth);
    points = plot(getIndexFns{fnn}(peakAlpha), expData{1,edn}, '*', 'MarkerSize', 8, 'LineWidth', lineWidth/2);
    bestIndexLabel = sprintf("%s_{−%s}", indexName{fnn}, as_4_dp_str(abs(peakAlpha)));
    pointsLabel = sprintf("%s and %s", bestIndexLabel, expData{2,edn});

    % Label the scatter plot
    title(sprintf('between %s and %s', expData{2,edn}, bestIndexLabel));
    xlabel(bestIndexLabel);
    ylabel(sprintf('%s', expData{2,edn}));
    xticklabels(strrep(xticklabels,'-','−'));
    yticklabels(strrep(yticklabels,'-','−'));
    leg = legend("Regression Line", "Actual Data");
    set(leg, 'location', "southeast");

    % Change the font size to size set in the start of the script
    set(findall(this_figure,'-property','FontSize'),'FontSize', fontSize)

    drawnow;
  end
end

if saveToFile
  % Save each figure to a separate file
  saveas(figure(1), "Scatter_E_T_F.png");
end
