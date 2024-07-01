% Combined_Corr_Curves_FTemp.m
% In this script, six values are closely approximated via golden section search
% - beta value for which correlation coefficient ρ is strongest between E_π and T_{b}
% E_π = Pi Electronic Energy
% T_{b} = F-Temeperature Index
% Additionally, curves ρ against α near these values are plotted in 1 figs

close all; % Close any figures already opened
clear;     % and clear all variables
format long; % More significant figures printed in the console
pkg load statistics;

% Config: Adjust lines' width, font size, whether or not to save to file
lineWidth = 2; % the line thickness in the graph
fontSize = 16; % the font size in the graph
saveToFile = true; % Set to true to auto-save plots

% Utility: Below is a function to round-off to 4 decimal places | returns string
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));

% Cell containing Entropy and Heat Capacity of 30 lower benzenoids
expData = {reshape([% Pi Electroni Energy
   8.0000  13.6832 19.3137 19.4483 24.9308 25.1875 25.1012 25.1922 25.2745 22.5055
   30.5440 30.7255 30.8805 30.8795 30.7627 30.9990 30.9362 30.9386 30.9432 30.8390
   30.9418 30.8336 28.2453 28.3361 28.2220 36.6814 31.4251 36.1557 34.5718 46.4974
   ]', 30, 1),
   "E_{π}" % Their labels
 };

% 30 by 3 array, number of edges in each dudv partition: (2,2), (2,3) then (3,3)
 d_f = [
   12 12 12 14 12 16 14 16 18 12 12 14 16 16 14 20 18 16 18 16 18 16 16 16 14 12 14 12 12 8
   0  4  8  6  12 8  10 8  6  8  16 14 12 12 14 8  10 12 10 12 10 10 8  8  10 12 10 20 12 16
   0  2  4  6  6  10 8  10 12 10 8  10 12 12 10 16 14 12 14 12 14 12 16 16 14 24 20 10 24 38
 ]'; % Used for computing indices based on edge endpoint degree partitions

% Define nu_exp here before using it in getIndexFns Lower 30 Bhs (nu_exp: number of vertices)
nu_exp = [
    6 10 14 14 18 18 18 18 18 18 22 22 22 22 22 22 22 22 22 22 22 22 20 20 20 24 22 26 24 32
]';

% Cell containing the three index-computing functions
% Usage: getIndexFns{n}(beta, nu_exp) | n=1:T_{b},
getIndexFns = {
    @(beta) d_f(:,1).*(2 ./ (nu_exp - 2)).^beta + d_f(:,2).*(2 ./ (nu_exp - 2) + 3 ./ (nu_exp - 3)).^beta + d_f(:,3).*(3 ./ (nu_exp - 3)).^beta;
}';

% Cell containing their respective labels
indexName = {"T"};

% Variables for indexing arrays and iterating for-loops
numData = size(expData, 2);       % one (E_{π})
numIndices = size(getIndexFns,1); % two (T_b)
numCurves = numData*numIndices;   % two (1 curves)

% All x in visible ranges (both plots - near and far)
xnear = linspace(-0.25, 0.2, 900);
xfar = linspace(-20,20,800); % for E_{π}

% Do the same procedure for each experimental data i.e., bp
for ii = 1:numData
  % This figure is for the zoomed-in plot
  figure(numData+ii); hold on;

  % Draw correlation curves for each index-property pair
  for n = 1:numIndices
    % Function to get corrcoef ρ between E_π with specified α
    %                                and T_{a}^1 /T_{a}^2 (depending on n)
    get_indices_vals = @(alpha) getIndexFns{n}(alpha)(!isnan(expData{1,ii}));
    ccFn = @(alpha) corrcoef(
      get_indices_vals(alpha), % Either T_{a}^1  or T_{a}^2
      expData{1,ii}(!isnan(expData{1,ii})) % E_{π}
    )(1,2);

    % generate corresponding y values
    ynear = arrayfun(ccFn, xnear);
    yfar = arrayfun(ccFn, xfar);

    % Compute peak values via. golden section search, and display in console
    disp(sprintf("%s against general_%s", expData{2,ii}, indexName{n}));
    peakAlpha = mean(GoldenSectionSearch_Maximum(ccFn, xnear(1), xnear(end), 1e-15));
    peakCorrCoeff = ccFn(peakAlpha);

    % Display peak alpha and peak correlation coefficient in the command window
    disp(['Peak Alpha: ', num2str(peakAlpha)]);
    disp(['Peak Correlation Coefficient: ', num2str(peakCorrCoeff)]);


    % Generate curve label [E_{π}] [T_{a}^1/T_{a}^2]
    curveLabels{n} = sprintf("%s and %s_B", expData{2,ii}, indexName{n});

    figure(ii); % One zoomed-out plot for each expData
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curvesFar(n) = plot(xfar, yfar, '-', 'LineWidth', lineWidth);
    drawnow;

    figure(numData+ii); % Each expData's set of curves has a zoomed-in plot
    hold on;
    % The curve is stored in a variable to be referenced in the legend spec
    curves(n) = plot(xnear, ynear, '-', 'LineWidth', lineWidth);

    % Show the peak in the plot: draw indicator lines & display coordinates
    plot([peakAlpha peakAlpha xnear(1)], [0 peakCorrCoeff peakCorrCoeff],
         '--k', 'LineWidth', lineWidth/2); % Black dashed indicator lines
    text(peakAlpha, peakCorrCoeff,
        {'', sprintf("(−%s, %s)", as_4_dp_str(abs(peakAlpha)), as_4_dp_str(peakCorrCoeff))},
        'VerticalAlignment', 'bottom');
        % Negative sign entered manually here to bypass the default
        % ... usage of "hypen-minus" instead of "minus" (− vs -)

    yend = min(ynear(1), ynear(end));
     ystart = max(ynear);% y value to be used as visible y lower bound
  end


  % Label this expData's zoomed-in plot
  xlabel('β');
  ylabel('ρ');
  leg = legend(curves, curveLabels); % curves contains all drawn "xnear" curves
  set(leg, 'location', "northeast"); % the location of the legend box

  ybox_space = 0.000005; % spacing between upper y-value with the blue dotted lines (figure 3)
  axis([xnear(1) xnear(end) yend ystart]); % Enforce figure's visible range
  drawnow;

  % Label the zoomed-out plot
  figure(ii);
  xlabel('β');
  ylabel('ρ');
  leg = legend(curvesFar, curveLabels); % curvesFar contains all drawn "xfar" curves
  set(leg, 'location', "southeast");

  hold off;
end

for ii = 1:2
  % Replace hyphens with minuses on negative axes
  figure(ii);
  xticklabels(strrep(xticklabels,'-','−'));
  yticklabels(strrep(yticklabels,'-','−'));

  % Set all fontsizes to size specified early in the script
  set(findall(figure(ii),'-property','FontSize'),'FontSize', fontSize)
end

if saveToFile
  saveas(figure(1), "01_comb_ccurves_E_{π}_FTemp_FAR.png");
  saveas(figure(2), "01_comb_ccurves_E_{π}_FTemp_NEAR.png");
end
