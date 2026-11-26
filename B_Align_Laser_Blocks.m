%% B_Align_Laser_Blocks
% Cuts the photometry, heart rate, RMSSD, pupil and running traces into
% laser pulse blocks. Time 0 is the onset (first pulse) of each block; each
% window runs from preBlockSec before the onset to postBlockSec after the
% end of the longest block.
%
% The two 465 nm traces are z-scored to their own baseline (baselineWinSec
% before time 0) in every block. Set normaliseTo405 = true to z-score the
% 405 nm traces the same way and subtract them (465z - 405z). Set
% zscorePupil = true to z-score the pupil to the same baseline.
%
% Any number of laser blocks is accepted (a warning is given if it is not
% expectedBlocks). If no laser pulses are found at all (e.g. a control
% day), the block times of the 25.11.2026 recording (laserTraceMat) are
% used instead, so the same periods of the recording are compared.
%
% Run after A_Parent_Align_Heka_Pupil_Photometry (uses the workspace), or
% on its own in a folder with one *_pupilAligned.mat file.
%
% Output (in the workspace, and saved to <HEKA name>_laserBlocks.mat)
%   blockTimeS     1 x nWin, time from block onset (s)
%   blockOnIdx     onset sample of each block (row of extractedMatCrop)
%   blockOffIdx    last laser-high sample of each block
%   blockDurS      duration of each block (s)
%   blocksFromTemplate  true if no laser was found and the template was used
%   blockMats      struct, one nBlocks x nWin matrix per signal:
%                  phot1_465z, phot2_465z, HR, RMSSD, pupil, run, laser
%                  (phot1_465z / phot2_465z hold 465z - 405z when
%                  normaliseTo405 = true; pupil is z-scored when
%                  zscorePupil = true)
%   blockPlotNames y-axis labels used in the figure

%% ---------------- Settings ----------------
normaliseTo405 = false;     % true = 465 z-score minus 405 z-score
zscorePupil    = false;     % true = z-score the pupil to the same baseline
preBlockSec    = 5;         % s before block onset
postBlockSec   = 5;         % s after block end
baselineWinSec = [-3 -1];   % s, z-score baseline relative to onset
blockGapSec    = 5;         % s, laser off for longer than this starts a new block
expectedBlocks = 6;         % warn if a different number of blocks is found
saveBlocks     = true;      % save the block matrices to a .mat
saveBlockFig   = true;      % save the figure as .fig and .png

% Used only if no laser pulses are found: block onsets and durations of
% 25112026__laserTraceMat.mat, in s from the first sample of extractedMat
% (the uncropped HEKA recording).
templateOnsetSec = [615.12 826.13 1037.1325 1248.13 1459.1325 1670.1325];
templateDurSec   = [28.49 28.49 28.49 28.4925 28.49 28.4925];

%% ---------------- Load the aligned data if not in the workspace ----------
if ~exist('extractedMatCrop', 'var')
    alignedFiles = dir('*_pupilAligned.mat');
    alignedFiles = alignedFiles(~startsWith({alignedFiles.name}, '._'));
    if numel(alignedFiles) ~= 1
        error(['Run A_Parent_Align_Heka_Pupil_Photometry first, or keep exactly ' ...
               'one *_pupilAligned.mat in:\n  %s'], pwd);
    end
    load(alignedFiles(1).name);
    fprintf('Loaded %s\n', alignedFiles(1).name);
end

% Title: mouse and day folder names, e.g. 'm947 d1 opto'
if ~exist('recordingTitle', 'var')
    [mousePath, dayName] = fileparts(pwd);
    [~, mouseName]       = fileparts(mousePath);
    recordingTitle       = [mouseName ' ' dayName];
end

nCrop     = size(extractedMatCrop, 1);
sampleSec = sampleIntervalMs / 1000;

%% ---------------- Find the laser blocks ----------------
laserCrop  = extractedMatCrop(:, 3);
laserEdges = diff(laserCrop);
pulseRise  = find(laserEdges == 1) + 1;     % first high sample of each pulse
pulseFall  = find(laserEdges == -1);        % last high sample of each pulse
if ~isempty(pulseRise)
    pulseFall = pulseFall(pulseFall >= pulseRise(1));
    if numel(pulseFall) < numel(pulseRise)  % laser still on at the end
        pulseFall(end + 1) = nCrop;
    end
    pulseFall = pulseFall(1:numel(pulseRise));
end

blocksFromTemplate = isempty(pulseRise);
if ~blocksFromTemplate
    pulseGapSec = (pulseRise(2:end) - pulseFall(1:end-1)) * sampleSec;
    newBlock    = [true; pulseGapSec > blockGapSec];
    lastOfBlock = [newBlock(2:end); true];
    blockOnIdx  = pulseRise(newBlock);
    blockOffIdx = pulseFall(lastOfBlock);
    fprintf('Laser: %d pulses in %d blocks\n', numel(pulseRise), numel(blockOnIdx));
else
    % No laser: use the template block times. They are in s from the first
    % sample of extractedMat; extractedMatCrop starts at row camRise.
    blockOnIdx  = round(templateOnsetSec(:) / sampleSec) + 1 - (camRise - 1);
    blockOffIdx = blockOnIdx + round(templateDurSec(:) / sampleSec) - 1;
    keepBlock   = blockOnIdx >= 1 & blockOnIdx <= nCrop;
    blockOnIdx  = blockOnIdx(keepBlock);
    blockOffIdx = blockOffIdx(keepBlock);
    warning(['No laser pulses found in %s: using the %d template block times ' ...
             'from 25112026__laserTraceMat.mat.'], recordingTitle, numel(blockOnIdx));
end
nBlocks   = numel(blockOnIdx);
blockDurS = (blockOffIdx - blockOnIdx + 1) * sampleSec;

if nBlocks == 0
    error('No laser blocks (or template blocks) fall within the recording of %s.', recordingTitle);
end
fprintf('%d blocks, duration %.2f to %.2f s\n', nBlocks, min(blockDurS), max(blockDurS));
if nBlocks ~= expectedBlocks
    warning('%s: %d laser blocks used, expected %d.', recordingTitle, nBlocks, expectedBlocks);
end

%% ---------------- Window around every block ----------------
preSamples  = round(preBlockSec / sampleSec);
postSamples = round((max(blockDurS) + postBlockSec) / sampleSec);
winOffsets  = -preSamples:postSamples;
blockTimeS  = winOffsets * sampleSec;
baseIdx     = blockTimeS >= baselineWinSec(1) & blockTimeS < baselineWinSec(2);

% Every signal is cut with the same window; samples outside the
% recording are NaN.
blockSignals = struct( ...
    'phot1_465', extractedMatCrop(:, 6), 'phot1_405', extractedMatCrop(:, 5), ...
    'phot2_465', extractedMatCrop(:, 9), 'phot2_405', extractedMatCrop(:, 8), ...
    'HR', hrCrop, 'RMSSD', hrvCrop.RMSSD, 'pupil', extractedMatCrop(:, 11), ...
    'run', runCrop, 'laser', laserCrop);
signalNames = fieldnames(blockSignals);

blockCut = struct();
for signalLoop = 1:numel(signalNames)
    fullTrace = blockSignals.(signalNames{signalLoop});
    cutMat    = nan(nBlocks, numel(winOffsets));
    for blockLoop = 1:nBlocks
        winIdx  = blockOnIdx(blockLoop) + winOffsets;
        inRange = winIdx >= 1 & winIdx <= nCrop;
        cutMat(blockLoop, inRange) = fullTrace(winIdx(inRange));
    end
    blockCut.(signalNames{signalLoop}) = cutMat;
end

%% ---------------- z-score to each block's own baseline ----------------
baseZ = @(cutMat) (cutMat - mean(cutMat(:, baseIdx), 2, 'omitnan')) ./ ...
                  std(cutMat(:, baseIdx), 0, 2, 'omitnan');

blockMats = struct();
if normaliseTo405
    blockMats.phot1_465z = baseZ(blockCut.phot1_465) - baseZ(blockCut.phot1_405);
    blockMats.phot2_465z = baseZ(blockCut.phot2_465) - baseZ(blockCut.phot2_405);
    photLabel = '465-405 (z)';
else
    blockMats.phot1_465z = baseZ(blockCut.phot1_465);
    blockMats.phot2_465z = baseZ(blockCut.phot2_465);
    photLabel = '465 (z)';
end
blockMats.HR    = blockCut.HR;
blockMats.RMSSD = blockCut.RMSSD;
if zscorePupil
    blockMats.pupil = baseZ(blockCut.pupil);
    pupilLabel = 'Pupil (z)';
else
    blockMats.pupil = blockCut.pupil;
    pupilLabel = 'Pupil (px)';
end
blockMats.run   = blockCut.run;
blockMats.laser = blockCut.laser;

%% ---------------- Plot ----------------
blockPlotFields = {'phot1_465z', 'phot2_465z', 'HR', 'RMSSD', 'pupil', 'run'};
blockPlotNames  = {['Phot1 ' photLabel], ['Phot2 ' photLabel], 'HR (bpm)', ...
                   'RMSSD (ms)', pupilLabel, 'Run (ticks/s)'};
blockFigTitle   = recordingTitle;
if blocksFromTemplate
    blockFigTitle = [recordingTitle ' (template block times)'];
end
blockFig = plotLaserBlockFigure(blockTimeS, blockMats, blockPlotFields, ...
                                blockPlotNames, max(blockDurS), blockFigTitle, false);

%% ---------------- Save ----------------
[~, hekaBase] = fileparts(hekaFile);
if saveBlocks
    blocksFile = fullfile(pwd, [hekaBase '_laserBlocks.mat']);
    save(blocksFile, 'blockMats', 'blockTimeS', 'blockOnIdx', 'blockOffIdx', ...
         'blockDurS', 'blocksFromTemplate', 'normaliseTo405', 'zscorePupil', ...
         'preBlockSec', 'postBlockSec', 'baselineWinSec', 'blockPlotFields', ...
         'blockPlotNames', 'recordingTitle', 'hekaFile');
    fprintf('Saved %s\n', blocksFile);
end
if saveBlockFig
    blockFigBase = fullfile(pwd, [hekaBase '_laserBlocks']);
    savefig(blockFig, [blockFigBase '.fig']);
    print(blockFig, [blockFigBase '.png'], '-dpng', '-r300');
    fprintf('Saved %s.fig and .png\n', blockFigBase);
end
