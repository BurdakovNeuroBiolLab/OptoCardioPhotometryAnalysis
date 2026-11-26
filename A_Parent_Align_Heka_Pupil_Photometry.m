%% A_Parent_Align_Heka_Pupil_Photometry
% Parent script: loads the pupil video data and the HEKA recording from the
% MATLAB current folder, puts the pupil diameter on the same time base as
% the photometry traces, and extracts heart rate and running.
%
% Same approach as A_Parent_Align_Heka_Pupil_Inscopix in
% OptoCardioInscopixAnalysis.
%
% Steps
%   1. pupilCircleFit.m          -> pupil diameter per video frame
%   2. A_Extract_HekaData.m      -> extractedMat (2.5 ms/sample, 10 columns)
%   3. alignPupilToPhotometry.m  -> pupil added as column 11, using the
%                                   Cam Synch high period (first rise to
%                                   last fall = start and end of the CSV)
%   4. Heart rate and HRV (RMSSD etc., DFA alpha1) from the ECG R peaks
%      (hrMethod: GRIN-code method or preprocessMultipleMice_2 method)
%   5. Running speed from the rotary encoder
%   6. Figure of photometry, laser, heart rate, RMSSD, pupil and running over the
%      camera period
%   7. Save
%
% The current folder must contain exactly one of each:
%   *zoompup*.csv (DeepLabCut pupil) and the HEKA export *L.mat.
%   Files starting with ._ are ignored.
%
% Output (in the workspace, and saved to <HEKA name>_pupilAligned.mat)
%   extractedMat      (nSamples x 11), 2.5 ms per row
%   extractedMatCrop  extractedMat cropped to the camera period
%   Columns: 1 HR (raw ECG)  5 Phot 1 405   9  Phot 2 465
%            2 Rotary 1      6 Phot 1 465   10 Cam synch
%            3 Laser         7 Phot 2 raw   11 Pupil diameter (px)
%            4 Phot 1 raw    8 Phot 2 405
%   colTitles         column names
%   hekaMs            time of every row of extractedMat (ms)
%   hekaMsCrop        time of every row of extractedMatCrop (ms)
%   pupilHekaT        centre of each video frame in HEKA time (ms)
%   hrCrop            heart rate (bpm) at every row of extractedMatCrop
%                     (beat-to-beat for hrMethod 'claude', windowed for 'nikMattia')
%   rT                R-peak times (HEKA s), refined to sub-sample precision
%   beatT, beatHR     time (HEKA s) and rate (bpm) of every beat-to-beat interval
%   okHR              true for beat intervals kept (not missed / extra beats)
%   hrvCrop           windowed HRV at every row of extractedMatCrop (struct):
%                     HR (bpm), RMSSD (ms), CVSD, CV_IBI, HRV_resid, nIBI
%   hrvMetrics        continuous_hrv_mouse output (window centres in HEKA s)
%   alpha1Whole       DFA alpha1 over all valid IBIs of the crop
%   runCrop           running (encoder ticks/s) at every row of extractedMatCrop

% C_Batch_Laser_Blocks sets batchMode = true and runs this script inside a
% function, so the workspace is already empty and must not be cleared.
if ~exist('batchMode', 'var') || ~batchMode
    clear; clc;
end

% Title for every figure: mouse and day folder names, e.g. 'm947 d1 opto'
[mousePath, dayName] = fileparts(pwd);
[~, mouseName]       = fileparts(mousePath);
recordingTitle       = [mouseName ' ' dayName];
figuresBefore        = findall(0, 'Type', 'figure');

%% ---------------- Settings ----------------
expectedFps = [];     % camera frame rate (Hz). If set, warn when the rate
                      % implied by the Cam Synch period differs by more
                      % than fpsTol. [] = no check.
fpsTol      = 0.02;   % fractional tolerance for the frame-rate check
saveAligned = true;   % save the aligned data to a .mat in the current folder
saveOverviewFig = true;  % save the final overview figure as .fig and .png

% Heart rate method
%   'claude'    : R peaks above a threshold, refined to sub-sample timing,
%                 missed / extra beats removed; HR is beat-to-beat
%                 (same as the GRIN code). Settings in the next block.
%   'nikMattia' : band-pass + findpeaks, no beat rejection; HR is the
%                 windowed HR from continuous_hrv_mouse (same as
%                 preprocessMultipleMice_2). Settings in the block after.
%                 Needs the Signal Processing Toolbox.
hrMethod = 'claude';

% 'claude' settings (defaults for mouse, as in the GRIN code)
hrRange       = [250 900]; % bpm; beats outside this range are dropped
hrThreshFrac  = 0.5;       % R-peak threshold as a fraction of the ECG 99.5th percentile
hrBaselineWin = 0.05;      % s, moving-mean window removed from the ECG
ibiRelTol     = 0.3;       % drop IBIs > 30% away from the running median (missed /
                           % extra beats)
showHRCheck   = false;     % true = extra figure with the detected R peaks
hrvWindowSec  = 10;        % s, sliding window for HRV (continuous_hrv_mouse)
hrvStepSec    = 0.25;      % s, step between HRV windows

% 'nikMattia' settings (as in preprocessMultipleMice_2)
nmEcgBand         = [8 30];  % Hz, band-pass
nmMinProm         = 0.15;    % findpeaks MinPeakProminence (ECG units - check
                             % this suits the HEKA ECG scale)
nmMaxWidth        = 0.15;    % s, findpeaks MaxPeakWidth
nmMinPeakDistance = 0.05;    % s, findpeaks MinPeakDistance (refractory period)
nmInvertECG       = false;   % true = flip the ECG before peak detection
nmHrvWindowSec    = 3;       % s, sliding window for HR / HRV
nmHrvStepSec      = 0.25;    % s, step between windows

% Running from the rotary encoder
runWinSec = 0.5;           % s, window over which encoder ticks are counted

%% ---------------- Put the child scripts on the path ----------------
alignRepoDir = fileparts(mfilename('fullpath'));
if ~isempty(alignRepoDir)
    addpath(alignRepoDir);
end
childScripts = {'pupilCircleFit', 'A_Extract_HekaData', 'alignPupilToPhotometry', ...
                'flag_bad_ibi', 'continuous_hrv_mouse', 'dfa_alpha1_local'};
for scriptLoop = 1:numel(childScripts)
    if isempty(which(childScripts{scriptLoop}))
        error(['Cannot find %s.m. Keep all the repository''s .m files in one ' ...
               'folder on the MATLAB path.'], childScripts{scriptLoop});
    end
end

pupilFiles = dir('*zoompup*.csv');
pupilFiles = pupilFiles(~startsWith({pupilFiles.name}, '._') & ...
                        ~endsWith({pupilFiles.name}, '_circlefit.csv'));
if isempty(pupilFiles)
    error('No pupil CSV (*zoompup*.csv) found in:\n  %s', pwd);
elseif numel(pupilFiles) > 1
    error('Found %d pupil CSVs (*zoompup*.csv) in:\n  %s\nKeep only one of:\n  %s', ...
          numel(pupilFiles), pwd, strjoin({pupilFiles.name}, sprintf('\n  ')));
end

hekaFiles = dir('*L.mat');
hekaFiles = hekaFiles(~startsWith({hekaFiles.name}, '._'));
if isempty(hekaFiles)
    error('No HEKA export (*L.mat) found in:\n  %s', pwd);
elseif numel(hekaFiles) > 1
    error('Found %d HEKA exports (*L.mat) in:\n  %s\nKeep only one of:\n  %s', ...
          numel(hekaFiles), pwd, strjoin({hekaFiles.name}, sprintf('\n  ')));
end

keepWorkspace = true;          % stops pupilCircleFit clearing our settings

%% ---------------- 1. Pupil ----------------
pupilCircleFit;                % defines diameter (cleaned, px), fname
pupilDiameter = diameter;
pupilFile     = fname;

%% ---------------- 2. HEKA ----------------
A_Extract_HekaData;            % defines extractedMat, colTitles, fileName
hekaFile = fileName;
fprintf('HEKA file: %s\n', hekaFile);

%% ---------------- 3. Pupil onto the HEKA time base ----------------
alignPupilToPhotometry;        % adds column 11, extractedMatCrop, camRise/Fall
hekaMsCrop = hekaMs(camRise:camFall);
hekaT      = hekaMsCrop / 1000;                    % crop samples, HEKA s
hekaFs     = 1000 / sampleIntervalMs;              % Hz
nCrop      = numel(hekaT);

if ~isempty(expectedFps)
    impliedFps = 1000 / pupilPeriodMs;
    if abs(impliedFps / expectedFps - 1) > fpsTol
        warning(['Cam Synch period implies %.2f Hz but expectedFps is %.2f Hz. ' ...
                 'Check that the pupil CSV and the HEKA file are from the same ' ...
                 'recording.'], impliedFps, expectedFps);
    end
end

%% ---------------- 4. Heart rate from the ECG ----------------
ecg         = extractedMatCrop(:, 1);
ibiRangeSec = sort(60 ./ hrRange);

switch hrMethod
case 'claude'
% Same method as A_Parent_Align_Heka_Pupil_Inscopix: R peaks are found in
% the full-rate ECG of the crop, refined to sub-sample precision and turned
% into beat-to-beat heart rate, which is interpolated onto every sample.
ecgHP    = ecg - movmean(ecg, max(3, round(hrBaselineWin * hekaFs)), 'omitnan');
if abs(localPrctile(ecgHP, 0.1)) > abs(localPrctile(ecgHP, 99.9))
    ecgHP = -ecgHP;                                        % R peaks point down
end
hrThresh = hrThreshFrac * localPrctile(ecgHP, 99.5);
candPeaks = find(ecgHP(2:end-1) > hrThresh & ecgHP(2:end-1) >= ecgHP(1:end-2) & ...
                 ecgHP(2:end-1) > ecgHP(3:end)) + 1;

% Refractory period: of two peaks closer than 60/maxBpm, keep the larger
minGap = round(60 / hrRange(2) * hekaFs);
rIdx = zeros(numel(candPeaks), 1); nR = 0;
for peakLoop = 1:numel(candPeaks)
    if nR > 0 && candPeaks(peakLoop) - rIdx(nR) < minGap
        if ecgHP(candPeaks(peakLoop)) > ecgHP(rIdx(nR))
            rIdx(nR) = candPeaks(peakLoop);
        end
    else
        nR = nR + 1; rIdx(nR) = candPeaks(peakLoop);
    end
end
rIdx = rIdx(1:nR);

% Parabolic interpolation through each peak and its two neighbours: the peak
% time is only known to +/- half a sample (1.25 ms at 400 Hz), which would
% otherwise add several bpm of jitter to beat-to-beat HR.
pkShift   = zeros(nR, 1);
innerPeak = rIdx > 1 & rIdx < numel(ecgHP);
ecgBefore = ecgHP(rIdx(innerPeak) - 1);
ecgPeak   = ecgHP(rIdx(innerPeak));
ecgAfter  = ecgHP(rIdx(innerPeak) + 1);
peakShift = 0.5 * (ecgBefore - ecgAfter) ./ (ecgBefore - 2 * ecgPeak + ecgAfter);
peakShift(~isfinite(peakShift)) = 0;
pkShift(innerPeak) = max(-0.5, min(0.5, peakShift));

rT     = hekaT(rIdx) + pkShift / hekaFs;                   % HEKA s
beatT  = (rT(1:end-1) + rT(2:end)) / 2;                    % HEKA s
beatHR = 60 ./ diff(rT);                                   % bpm
okHR   = ~flag_bad_ibi(diff(rT), ibiRangeSec, ibiRelTol, 21);

% HRV: missed / extra beats excluded, as above
hrvWinUsed   = hrvWindowSec;
hrvStepUsed  = hrvStepSec;
hrvRejectArg = {ibiRangeSec, ibiRelTol};

case 'nikMattia'
% Same method as preprocessMultipleMice_2: band-pass the ECG, find R peaks
% with findpeaks, no beat rejection. HR is the windowed HR from
% continuous_hrv_mouse (3 s windows by default).
if nmInvertECG
    ecg = -ecg;
end
ecgHP = bandpass(ecg, nmEcgBand, hekaFs);
[~, peakLocsSec] = findpeaks(ecgHP, hekaFs, 'MaxPeakWidth', nmMaxWidth, ...
    'MinPeakProminence', nmMinProm, 'MinPeakDistance', nmMinPeakDistance);
peakLocsSec = peakLocsSec(:);
rIdx     = round(peakLocsSec * hekaFs) + 1;                % sample in the crop
nR       = numel(rIdx);
hrThresh = NaN;                                            % no fixed threshold
rT       = hekaT(1) + peakLocsSec;                         % HEKA s
beatT    = (rT(1:end-1) + rT(2:end)) / 2;                  % HEKA s
beatHR   = 60 ./ diff(rT);                                 % bpm
okHR     = true(size(beatHR));                             % no beat rejection

% HRV: no beat rejection, as in preprocessMultipleMice_2
hrvWinUsed   = nmHrvWindowSec;
hrvStepUsed  = nmHrvStepSec;
hrvRejectArg = {};

otherwise
    error('Unknown hrMethod ''%s''. Use ''claude'' or ''nikMattia''.', hrMethod);
end

%% ---------------- 4b. HRV ----------------
% Windowed HRV from continuous_hrv_mouse, interpolated onto every sample
% of the crop.
hrvFields  = {'HR', 'RMSSD', 'CVSD', 'CV_IBI', 'HRV_resid', 'nIBI'};
hrvCrop    = struct();
for fieldLoop = 1:numel(hrvFields)
    hrvCrop.(hrvFields{fieldLoop}) = nan(nCrop, 1);
end
hrvMetrics = struct();
if nR >= 3 && rT(end) - rT(1) > hrvWinUsed
    hrvMetrics = continuous_hrv_mouse(rT, hrvWinUsed, hrvStepUsed, hrvRejectArg{:});
    for fieldLoop = 1:numel(hrvFields)
        hrvCrop.(hrvFields{fieldLoop}) = interp1(hrvMetrics.t, ...
            hrvMetrics.(hrvFields{fieldLoop}), hekaT, 'linear', NaN);
    end
    hrvCrop.RMSSD = 1000 * hrvCrop.RMSSD;                   % s -> ms
else
    warning('Recording shorter than the HRV window (%g s); HRV not computed.', hrvWinUsed);
end

%% ---------------- 4c. Heart rate on every sample ----------------
% 'claude': beat-to-beat HR; 'nikMattia': windowed HR
if strcmp(hrMethod, 'claude')
    if nnz(okHR) >= 2
        hrCrop = interp1(beatT(okHR), beatHR(okHR), hekaT, 'linear', NaN);
    else
        hrCrop = nan(nCrop, 1);
        warning('Fewer than 2 heart beats detected; check hrThreshFrac / hrRange.');
    end
else
    hrCrop = hrvCrop.HR;
end
fprintf('Heart rate (%s): %d R peaks, median %.0f bpm (%d IBIs dropped)\n', ...
        hrMethod, nR, median(beatHR(okHR)), nnz(~okHR));

% DFA alpha1 once over the whole crop (too few beats per window to be
% reliable). Uses the valid IBIs only, joined across rejected ones.
ibiValid    = diff(rT);
ibiValid    = ibiValid(okHR);
alpha1Whole = dfa_alpha1_local(ibiValid);
fprintf('HRV: %g s windows, median RMSSD %.2f ms; DFA alpha1 (whole crop) %.2f\n', ...
        hrvWinUsed, median(hrvCrop.RMSSD, 'omitnan'), alpha1Whole);

%% ---------------- 5. Running from the rotary encoder ----------------
% Every state change of the binarised encoder is one tick; ticks are
% counted in a sliding window of runWinSec and converted to ticks per second.
% The window is trailing (the runWinSec before each sample), so running
% never appears before it happened, e.g. before a laser onset.
runTicks = [0; abs(diff(extractedMatCrop(:, 2)))];
runWinSamples = max(1, round(runWinSec * hekaFs));
runCrop  = movsum(runTicks, [runWinSamples-1 0]) / (runWinSamples / hekaFs);

%% ---------------- 6. Final figure over the camera period ----------------
cropTimeS = hekaT - hekaT(1);                        % s from camera start

plotTraces = {extractedMatCrop(:, 6), extractedMatCrop(:, 5), ...
              extractedMatCrop(:, 9), extractedMatCrop(:, 8), ...
              extractedMatCrop(:, 3), hrCrop, hrvCrop.RMSSD, ...
              extractedMatCrop(:, 11), runCrop};
plotNames  = {'Phot1 465', 'Phot1 405', 'Phot2 465', 'Phot2 405', ...
              'Laser', 'HR (bpm)', 'RMSSD (ms)', 'Pupil (px)', 'Run (ticks/s)'};
plotColour = {[0 0.6 0], [0.7 0 0.7], [0 0.6 0], [0.7 0 0.7], ...
              [0 0 0.8], [0.85 0.1 0.1], [0.9 0.4 0.7], [0.9 0.5 0], [0.2 0.6 0.9]};

nPlots = numel(plotTraces);
overviewFig = figure('Name', 'Photometry, heart rate, pupil and running', ...
                     'Renderer', 'painters', 'Position', [900 50 700 900]);
axPlots = gobjects(nPlots, 1);
for plotLoop = 1:nPlots
    axPlots(plotLoop) = subplot(nPlots, 1, plotLoop);
    plot(cropTimeS, plotTraces{plotLoop}, 'Color', plotColour{plotLoop})
    ylabel(plotNames{plotLoop})
    box off
    if plotLoop < nPlots
        set(axPlots(plotLoop), 'XTickLabel', []);
    end
end
xlabel(axPlots(end), 'Time from camera start (s)')
linkaxes(axPlots, 'x');
xlim(axPlots(end), [cropTimeS(1) cropTimeS(end)]);

%% ---------------- Optional: ECG peak check ----------------
if showHRCheck
    figure('Name', 'ECG R-peak detection', 'Position', [100 100 1200 400]);
    plot(cropTimeS, ecgHP, 'k'); hold on
    plot(rT - hekaT(1), ecgHP(rIdx), 'r.', 'MarkerSize', 10);
    if ~isnan(hrThresh), yline(hrThresh, 'b--'); end
    hold off
    xlabel('Time from camera start (s)'); ylabel('ECG (baseline removed)');
end

%% ---------------- Title every figure made by this script ----------------
newFigures = setdiff(findall(0, 'Type', 'figure'), figuresBefore);
for titleLoop = 1:numel(newFigures)
    sgtitle(newFigures(titleLoop), recordingTitle, 'Interpreter', 'none');
end

%% ---------------- 7. Save ----------------
[~, hekaBase] = fileparts(hekaFile);
if saveOverviewFig
    overviewBase = fullfile(pwd, [hekaBase '_overview']);
    savefig(overviewFig, [overviewBase '.fig']);
    print(overviewFig, [overviewBase '.png'], '-dpng', '-r300');
    fprintf('Saved %s.fig and .png\n', overviewBase);
end
if saveAligned
    alignedFile = fullfile(pwd, [hekaBase '_pupilAligned.mat']);
    save(alignedFile, 'extractedMat', 'extractedMatCrop', 'colTitles', ...
         'hekaMs', 'hekaMsCrop', 'camRise', 'camFall', 'pupilDiameter', ...
         'pupilHekaT', 'pupilPeriodMs', 'sampleIntervalMs', ...
         'hrCrop', 'rT', 'beatT', 'beatHR', 'okHR', 'runCrop', 'runWinSec', ...
         'hrvCrop', 'hrvMetrics', 'alpha1Whole', 'hrMethod', 'hrvWinUsed', 'hrvStepUsed', ...
         'hekaFile', 'pupilFile', 'recordingTitle');
    fprintf('Saved %s\n', alignedFile);
end

%% ---------------- Helper ----------------
function prctileValue = localPrctile(values, pct)
% Percentile without the Statistics Toolbox
values = sort(values(~isnan(values(:))));
if isempty(values), prctileValue = 0; return; end
prctileValue = values(max(1, min(numel(values), round(pct / 100 * numel(values)))));
end
