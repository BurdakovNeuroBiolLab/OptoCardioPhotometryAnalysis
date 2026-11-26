function figHandle = plotLaserBlockFigure(timeS, traceMats, plotFields, plotNames, ...
                                          laserDurS, figTitle, showMedian)
% PLOTLASERBLOCKFIGURE  Laser-block figure used by B_Align_Laser_Blocks and
% C_Batch_Laser_Blocks.
%
%   timeS       1 x nTime, time from block onset (s)
%   traceMats   struct with one nTraces x nTime matrix per field
%   plotFields  fields of traceMats to plot, one subplot each (top to bottom)
%   plotNames   y-axis labels, same order as plotFields
%   laserDurS   laser period (s) shaded in light blue from time 0
%   figTitle    title above the figure (e.g. 'm947 d1 opto')
%   showMedian  true = also plot the median across traces (thick dotted black)
%
% Every row of each matrix is plotted in its own colour; the mean across
% rows is a thick black line.

if nargin < 7
    showMedian = false;
end

nPlots    = numel(plotFields);
nTraces   = size(traceMats.(plotFields{1}), 1);
traceCols = lines(max(nTraces, 1));

figHandle = figure('Renderer', 'painters', 'Position', [500 200 400 800]);
axPlots   = gobjects(nPlots, 1);
for plotLoop = 1:nPlots
    axPlots(plotLoop) = subplot(nPlots, 1, plotLoop);
    hold on
    cutMat = traceMats.(plotFields{plotLoop});

    % Laser period shaded in light blue
    yRange = [min(cutMat(:), [], 'omitnan'), max(cutMat(:), [], 'omitnan')];
    if isempty(yRange) || any(isnan(yRange)) || yRange(1) == yRange(2)
        yRange = [0 1];
    end
    patch([0 laserDurS laserDurS 0], yRange([1 1 2 2]), ...
          [0.85 0.9 1], 'EdgeColor', 'none', 'HandleVisibility', 'off');

    for traceLoop = 1:nTraces
        plot(timeS, cutMat(traceLoop, :), 'Color', traceCols(traceLoop, :), 'LineWidth', 0.5);
    end
    plot(timeS, mean(cutMat, 1, 'omitnan'), 'k', 'LineWidth', 2);
    if showMedian
        plot(timeS, median(cutMat, 1, 'omitnan'), 'k:', 'LineWidth', 2);
    end
    hold off

    ylabel(plotNames{plotLoop})
    box off
    if plotLoop < nPlots
        set(axPlots(plotLoop), 'XTickLabel', []);
    end
end
xlabel(axPlots(end), 'Time from block onset (s)')
linkaxes(axPlots, 'x');
xlim(axPlots(end), [timeS(1) timeS(end)]);
sgtitle(figHandle, figTitle, 'Interpreter', 'none');
end
