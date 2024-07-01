% Good_alpha_intervals_FTemp.m
% This script generates six plots showing beta intervals for:
% - Good correlation coefficient ρ between E_{π} and T_F

close all; % Close any open figures
clear;     % Clear all workspace variables
format long; % Set format for long floating-point numbers

% Configurations (Adjust line width, font size, and save-to-file option)
lineWidth = 3;
fontSize = 24;
saveToFile = false; % Set to true to save plots automatically

% Utility function to round to 4 decimal places and convert to string
as_4_dp_str = @(x) sprintf('%.04f', round(x*(10^4))/(10^4));

% Experimental data: Entropy and Heat Capacity of 22 lower benzenoids
expData = {reshape([% Pi Electronic Energy
   8.0000  13.6832 19.3137 19.4483 24.9308 25.1875 25.1012 25.1922 25.2745 22.5055
   30.5440 30.7255 30.8805 30.8795 30.7627 30.9990 30.9362 30.9386 30.9432 30.8390
   30.9418 30.8336 28.2453 28.3361 28.2220 36.6814 31.4251 36.1557 34.5718 46.4974
   ]', 30, 1),
   "E_{π}" % Label
 };

% Number of edges in each dudv partition: (2,2), (2,3), then (3,3)
d_f = [
   12 12 12 14 12 16 14 16 18 12 12 14 16 16 14 20 18 16 18 16 18 16 16 16 14 12 14 12 12 8
   0  4  8  6  12 8  10 8  6  8  16 14 12 12 14 8  10 12 10 12 10 10 8  8  10 12 10 20 12 16
   0  2  4  6  6  10 8  10 12 10 8  10 12 12 10 16 14 12 14 12 14 12 16 16 14 24 20 10 24 38
 ]';

% Number of vertices (nu_exp) for lower benzenoids
nu_exp = [
    6 10 14 14 18 18 18 18 18 18 22 22 22 22 22 22 22 22 22 22 22 22 20 20 20 24 22 26 24 32
]';

% Index-computing function for T_F
getIndexFns = {
    @(beta) d_f(:,1) .* (2 ./ (nu_exp - 2)).^beta + ...
             d_f(:,2) .* (2 ./ (nu_exp - 2) + 3 ./ (nu_exp - 3)).^beta + ...
             d_f(:,3) .* (3 ./ (nu_exp - 3)).^beta;
}';

% Index labels
indexName = {"T_B"};

% Number of experimental data and indices
numData = size(expData, 2);
numIndices = size(getIndexFns, 1);
numCurves = numData * numIndices;

% Bounds for visible intervals for each index-property pair
xstart = [-2]; % Start of interval for E_{π}
xend = [0.7];  % End of interval for E_{π}

% Exact ρ value considered good for E_{π}
a_goodrho = [0.98];

% Colors for shading and indicators
colShaded = {[0.85, 1, 1]};
colIndicatorVert = {[0.2, 0.55, 0.55]};
colIndicatorHorz = {[0.35, 0.75, 0.75]};
colCurve = {[0, 0.75, 0.75]};

% Iterate over each experimental data (E_{π})
for ii = 1:numData
    % Create a new figure for each experimental data
    figure('Name', sprintf('Correlation between %s and %s', indexName{1}, expData{1,ii}));
    hold on;

    % Plot correlation curves for each index-property pair
    for n = 1:numIndices
        % Function to compute correlation coefficient ρ between E_{π} and T_F
        ccFn = @(beta) corrcoef(...
            getIndexFns{n}(beta), ...
            expData{1,ii} ...
        )(1,2);

        % Find beta for maximum ρ
        peakBeta = mean(GoldenSectionSearch_Maximum(ccFn, -4, 0, 1e-15));

        % Function to find interval limits where ρ is close to a_goodrho
        ccFn_good = @(b) -abs(ccFn(b) - a_goodrho(ii));
        getLimitFromInterval = @(lb, ub) mean(GoldenSectionSearch_Maximum(ccFn_good, lb, ub, 1e-15));

        % Compute interval limits
        b_lb = getLimitFromInterval(peakBeta - 3, peakBeta);
        b_ub = getLimitFromInterval(peakBeta, peakBeta + 3);

        % Display interval limits in console
        disp(sprintf("ρ(%s, %s) ≥ %.02f when β ∈ [%.08f, %.08f]", expData{2,ii}, indexName{n}, a_goodrho(ii), b_lb, b_ub));

        % Plot the curve excluding the good beta range (to be separately drawn)
        x = [linspace(xstart(ii), b_lb, 300), linspace(b_ub, xend(ii), 300)];
        y = arrayfun(ccFn, x);
        plot(x, y, '-', 'LineWidth', lineWidth);

        % Shade the area inside the good beta interval
        rectangle('Position', [b_lb, min(y), b_ub - b_lb, max(y) - min(y)], ...
                  'FaceColor', colShaded{ii}, 'LineStyle', 'none');

        % Draw indicator lines for the good beta interval
        plot([b_lb, b_lb], [min(y), max(y)], '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});
        plot([b_ub, b_ub], [min(y), max(y)], '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorVert{ii});

        % Horizontal dashed lines
        plot([xstart(ii), b_lb], [a_goodrho(ii), a_goodrho(ii)], '--k', 'LineWidth', lineWidth/1.75);
        plot([b_lb, b_ub], [a_goodrho(ii), a_goodrho(ii)], '--', 'LineWidth', lineWidth/1.75, 'Color', colIndicatorHorz{ii});

        % Label the plot
        title(sprintf('Correlation between %s and %s', expData{2,ii}, indexName{n}));
        xlabel('β');
        ylabel('ρ');
        drawnow;
        axis([xstart(ii), xend(ii), min(y), max(y)]);

        % Replace all hyphens with minus signs
        xticklabels(strrep(xticklabels,'-','−'));
        yticklabels(strrep(yticklabels,'-','−'));

        % Set font size
        set(gca, 'FontSize', fontSize);
        drawnow;

        hold off;
    end
end

% Save figures if saveToFile is true
if saveToFile
    saveas(figure(1), 'correlation_T_F_E_pi.png');
end
