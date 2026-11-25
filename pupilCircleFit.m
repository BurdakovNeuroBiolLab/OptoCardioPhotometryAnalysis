% pupilCircleFit.m
% Finds a DeepLabCut CSV with "zoompup" in its name in the current folder,
% removes low-confidence points (set to NaN), fits a circle to the remaining
% pupil points in every frame, removes blinks / bad fits / outliers, and
% plots radius and diameter over time (before vs after cleaning).
%
% Same logic as pupilCircleFit.m in OptoCardioInscopixAnalysis.
%
% File layout: 3 header lines, then per row:
%   col 1 = frame index / time, then repeating [x, y, likelihood] per point.
%
% Cleaning steps:
%   1. Fit only frames with >= minPoints good points
%   2. Remove blinks (frames with many low-confidence points) +/- blinkPad frames
%   3. Remove poor fits (point scatter around the circle > maxRelResid * radius)
%   4. (optional) Remove remaining single-frame spikes with a moving-median test
%   5. (optional) Fill short NaN gaps by linear interpolation
%
% Outputs used by the parent script
%   diameter   cleaned pupil diameter per video frame (pixels)
%   fname      full path of the pupil CSV

% When called from a parent script that sets keepWorkspace = true, the
% workspace is not cleared so the parent's variables survive.
if ~exist('keepWorkspace', 'var') || ~keepWorkspace
    clear; clc;
end

%% Settings
confThresh   = 0.99;  % points with likelihood < confThresh are set to NaN
minPoints    = 6;     % step 1: minimum good points needed to fit a circle
blinkMinBad  = 3;     % step 2: frame is a blink if >= this many points are low-confidence
blinkPad     = 2;     % step 2: also remove this many frames before/after each blink
maxRelResid  = 0.1;   % step 3: reject fits where RMS point scatter > this fraction of radius
doOutliers   = true;  % step 4: remove spikes that deviate from the moving median
outlierWin   = 31;    % step 4: moving-median window (frames)
outlierK     = 5;     % step 4: flag if > outlierK scaled MADs from the moving median
doFillGaps   = false; % step 5: interpolate across short NaN gaps
maxGapFill   = 10;    % step 5: longest gap (frames) to fill
fps          = [];    % camera frame rate (Hz). Leave empty to plot vs frame number.

%% Find and load the file
pupilCsvFiles = dir(fullfile(pwd, '*zoompup*.csv'));
% Ignore macOS "._" ghost files and this script's own output files
pupilCsvNames = {pupilCsvFiles.name};
keepFile      = ~startsWith(pupilCsvNames, '._') & ~endsWith(pupilCsvNames, '_circlefit.csv');
pupilCsvFiles = pupilCsvFiles(keepFile);
if isempty(pupilCsvFiles)
    error('No CSV file containing "zoompup" found in %s', pwd);
end
if numel(pupilCsvFiles) > 1
    warning('%d "zoompup" files found; using the first: %s', numel(pupilCsvFiles), pupilCsvFiles(1).name);
end
fname = fullfile(pupilCsvFiles(1).folder, pupilCsvFiles(1).name);
fprintf('Loading %s\n', pupilCsvFiles(1).name);

csvData = readmatrix(fname, 'NumHeaderLines', 3);

frame   = csvData(:, 1);
nFrames = size(csvData, 1);
nPoints = floor((size(csvData, 2) - 1) / 3);

xCols     = 2:3:(1 + 3*nPoints);
pointX    = csvData(:, xCols);          % frames x points
pointY    = csvData(:, xCols + 1);
pointConf = csvData(:, xCols + 2);

%% Replace low-confidence points with NaN
badPoint = pointConf < confThresh;
pointX(badPoint) = NaN;
pointY(badPoint) = NaN;
fprintf('Points set to NaN (likelihood < %g): %.1f%%\n', confThresh, 100*mean(badPoint(:)));

%% Fit a circle to each frame (any frame with >= 3 good points)
xc       = nan(nFrames, 1);
yc       = nan(nFrames, 1);
radius   = nan(nFrames, 1);
relResid = nan(nFrames, 1);
nValid   = sum(~badPoint, 2);

for frameLoop = 1:nFrames
    goodPoint = ~isnan(pointX(frameLoop, :)) & ~isnan(pointY(frameLoop, :));
    if nnz(goodPoint) >= 3
        framePointsX = pointX(frameLoop, goodPoint)';
        framePointsY = pointY(frameLoop, goodPoint)';
        [xc(frameLoop), yc(frameLoop), radius(frameLoop)] = fitCircle(framePointsX, framePointsY);
        relResid(frameLoop) = sqrt(mean((hypot(framePointsX - xc(frameLoop), ...
            framePointsY - yc(frameLoop)) - radius(frameLoop)).^2)) / radius(frameLoop);
    end
end
radiusRaw = radius;             % uncleaned trace, plotted in black

%% Step 1: require enough good points
rejFewPts = nValid < minPoints & ~isnan(radius);
radius(rejFewPts) = NaN;
fprintf('Step 1: frames removed with < %d good points: %d\n', minPoints, nnz(rejFewPts));

%% Step 2: blinks (+ padding)
blink = sum(badPoint, 2) >= blinkMinBad;
blink = movmax(double(blink), 2*blinkPad + 1) > 0;
rejBlink = blink & ~isnan(radius);
radius(blink) = NaN;
fprintf('Step 2: frames removed as blinks (+/- %d frames): %d\n', blinkPad, nnz(rejBlink));

%% Step 3: poor circle fits
rejResid = relResid > maxRelResid & ~isnan(radius);
radius(rejResid) = NaN;
fprintf('Step 3: frames removed for poor fit (rel. residual > %g): %d\n', ...
        maxRelResid, nnz(rejResid));

%% Step 4: remaining spikes (moving-median outlier test)
rejOutlier = false(nFrames, 1);
if doOutliers
    rejOutlier = isoutlier(radius, 'movmedian', outlierWin, ...
                           'ThresholdFactor', outlierK) & ~isnan(radius);
    radius(rejOutlier) = NaN;
end
fprintf('Step 4: frames removed as moving-median outliers: %d\n', nnz(rejOutlier));

%% Step 5: fill short gaps
if doFillGaps
    nNaNBefore = nnz(isnan(radius));
    radius     = fillmissing(radius, 'linear', 'MaxGap', maxGapFill);
    fprintf('Step 5: NaN frames filled (gaps <= %d frames): %d\n', ...
            maxGapFill, nNaNBefore - nnz(isnan(radius)));
end

% Keep the centre consistent with the cleaned radius
xc(isnan(radius) & ~isnan(radiusRaw)) = NaN;
yc(isnan(radius) & ~isnan(radiusRaw)) = NaN;

diameterRaw = 2 * radiusRaw;
diameter    = 2 * radius;
fprintf('Total frames without a valid radius: %d of %d (%.2f%%)\n', ...
        nnz(isnan(radius)), nFrames, 100*mean(isnan(radius)));

%% Time axis
if isempty(fps)
    frameTime = frame;       frameTimeLabel = 'Frame';
else
    frameTime = frame / fps; frameTimeLabel = 'Time (s)';
end

%% Figure 1: radius and diameter, before (black) vs after cleaning (colour)
figure('Name', 'Pupil size over time');
axRadius = subplot(2, 1, 1);
plot(frameTime, radiusRaw, 'k-'); hold on;
plot(frameTime, radius, 'b-');    hold off;
ylabel('Radius (pixels)');
title(sprintf('Pupil circle fit  (likelihood < %g set to NaN)', confThresh));
legend({'before cleaning', 'after cleaning'}, 'Location', 'best');
grid on;

axDiameter = subplot(2, 1, 2);
plot(frameTime, diameterRaw, 'k-'); hold on;
plot(frameTime, diameter, 'r-');    hold off;
ylabel('Diameter (pixels)');
xlabel(frameTimeLabel);
legend({'before cleaning', 'after cleaning'}, 'Location', 'best');
grid on;
linkaxes([axRadius axDiameter], 'x');

%% Figure 2: check the fit on the first valid frame
checkFrame = find(~isnan(radius), 1);
theta = linspace(0, 2*pi, 200);
figure('Name', 'Circle fit check');
hold on;
plot(csvData(checkFrame, xCols), csvData(checkFrame, xCols + 1), 'o', 'Color', [0.7 0.7 0.7]); % all points
plot(pointX(checkFrame, :), pointY(checkFrame, :), 'ko', 'MarkerFaceColor', 'k');            % points used
plot(xc(checkFrame) + radius(checkFrame)*cos(theta), ...
     yc(checkFrame) + radius(checkFrame)*sin(theta), 'r-', 'LineWidth', 1.5);
plot(xc(checkFrame), yc(checkFrame), 'r+');
hold off;
axis equal; set(gca, 'YDir', 'reverse'); grid on; box on;
xlabel('x (pixels)'); ylabel('y (pixels)');
legend({'all points', 'used in fit', 'fitted circle', 'centre'}, 'Location', 'best');
title(sprintf('Frame %g: r = %.2f px, d = %.2f px', frame(checkFrame), ...
      radius(checkFrame), diameter(checkFrame)));

%% Results table (saved alongside the CSV)
results = table(frame, xc, yc, radiusRaw, diameterRaw, radius, diameter, nValid, ...
    relResid, rejFewPts, rejBlink, rejResid, rejOutlier, ...
    'VariableNames', {'frame', 'centerX', 'centerY', 'radiusRaw', 'diameterRaw', ...
    'radius', 'diameter', 'nValidPoints', 'relResidual', ...
    'removedFewPoints', 'removedBlink','removedPoorFit', 'removedOutlier'});
[~, pupilCsvBase] = fileparts(fname);
outName = fullfile(pupilCsvFiles(1).folder, [pupilCsvBase '_circlefit.csv']);
writetable(results, outName);
fprintf('Saved results to %s\n', outName);

%% ---------------------------------------------------------------------
function [centreX, centreY, radiusFit] = fitCircle(pointsX, pointsY)
% Algebraic least-squares circle fit (Kasa method).
% Solves x^2 + y^2 + D*x + E*y + F = 0 for D, E, F.
    designMat = [pointsX, pointsY, ones(size(pointsX))];
    rhsVec    = -(pointsX.^2 + pointsY.^2);
    coeffs    = designMat \ rhsVec;
    centreX   = -coeffs(1) / 2;
    centreY   = -coeffs(2) / 2;
    radiusFit = sqrt(centreX^2 + centreY^2 - coeffs(3));
end
