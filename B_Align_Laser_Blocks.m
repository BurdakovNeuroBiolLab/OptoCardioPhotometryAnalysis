%% B_Align_Laser_Blocks
% Cuts the photometry, heart rate, RMSSD, pupil and running traces into
% laser pulse blocks. Time 0 is the onset (first pulse) of each block; each
% window runs from preBlockSec before the onset to postBlockSec after the
% end of the longest block.
%
% The two 465 nm traces are z-scored to their own baseline (baselineWinSec
% before time 0) in every block. Set normaliseTo405 = true to z-score the
% 405 nm traces the same way and subtract them (465z - 405z).
%
% Run after A_Parent_Align_Heka_Pupil_Photometry (uses the workspace), or
% on its own in a folder with one *_pupilAligned.mat file.
%
% Output (in the workspace, and saved to <HEKA name>_laserBlocks.mat)
%   blockTimeS     1 x nWin, time from block onset (s)
%   blockOnIdx     onset sample of each block (row of extractedMatCrop)
%   blockOffIdx    last laser-high sample of each block
%   blockDurS      duration of each block (s)
%   blockMats      struct, one nBlocks x nWin matrix per signal:
%                  phot1_465z, phot2_465z, HR, RMSSD, pupil, run, laser
%                  (phot1_465z / phot2_465z hold 465z - 405z when
%                  normaliseTo405 = true)

%% ---------------- Settings ----------------
normaliseTo405 = false;     % true = 465 z-score minus 405 z-score
preBlockSec    = 5;         % s before block onset
postBlockSec   = 5;         % s after block end
baselineWinSec = [-3 -1];   % s, z-score baseline relative to onset
blockGapSec    = 5;         % s, laser off for longer than this starts a new block
expectedBlocks = 6;         % warn if a different number of blocks is found
saveBlocks     = true;      % save the block matrices to a .mat

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

nCrop     = size(extractedMatCrop, 1);
sampleSec = sampleIntervalMs / 1000;

%% ---------------- Find the laser blocks ----------------
laserCrop   = extractedMatCrop(:, 3);
laserEdges  = diff(laserCrop);
pulseRise   = find(laserEdges == 1) + 1;     % first high sample of each pulse
pulseFall   = find(laserEdges == -1);        % last high sample of each pulse
pulseFall   = pulseFall(pulseFall >= pulseRise(1));
pulseRise   = pulseRise(1:numel(pulseFall));

pulseGapSec = (pulseRise(2:end) - pulseFall(1:end-1)) * sampleSec;
newBlock    = [true; pulseGapSec > blockGapSec];
lastOfBlock = [newBlock(2:end); true];
blockOnIdx  = pulseRise(newBlock);
blockOffIdx = pulseFall(lastOfBlock);
nBlocks     = numel(blockOnIdx);
blockDurS   = (blockOffIdx - blockOnIdx + 1) * sampleSec;

fprintf('Laser: %d pulses in %d blocks (block duration %.2f to %.2f s)\n', ...
        numel(pulseRise), nBlocks, min(blockDurS), max(blockDurS));
if nBlocks ~= expectedBlocks
    warning('Found %d laser blocks, expected %d. Check blockGapSec.', nBlocks, expectedBlocks);
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

%% ---------------- z-score photometry to each block's own baseline -------
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
blockMats.pupil = blockCut.pupil;
blockMats.run   = blockCut.run;
blockMats.laser = blockCut.laser;

%% ---------------- Plot ----------------
plotFields = {'phot1_465z', 'phot2_465z', 'HR', 'RMSSD', 'pupil', 'run'};
plotNames  = {['Phot1 ' photLabel], ['Phot2 ' photLabel], 'HR (bpm)', ...
              'RMSSD (ms)', 'Pupil (px)', 'Run (ticks/s)'};
blockColours = lines(nBlocks);
nPlots       = numel(plotFields);

figure('Renderer', 'painters', 'Position', [500 200 400 800]);
axPlots = gobjects(nPlots, 1);
for plotLoop = 1:nPlots
    axPlots(plotLoop) = subplot(nPlots, 1, plotLoop);
    hold on
    cutMat = blockMats.(plotFields{plotLoop});

    % Laser period (shortest to longest block) shaded in light blue
    yRange = [min(cutMat(:), [], 'omitnan'), max(cutMat(:), [], 'omitnan')];
    if any(isnan(yRange)) || yRange(1) == yRange(2)
        yRange = [0 1];
    end
    patch([0 max(blockDurS) max(blockDurS) 0], yRange([1 1 2 2]), ...
          [0.85 0.9 1], 'EdgeColor', 'none', 'HandleVisibility', 'off');

    for blockLoop = 1:nBlocks
        plot(blockTimeS, cutMat(blockLoop, :), ...
             'Color', blockColours(blockLoop, :), 'LineWidth', 0.5);
    end
    plot(blockTimeS, mean(cutMat, 1, 'omitnan'), 'k', 'LineWidth', 2);
    hold off

    ylabel(plotNames{plotLoop})
    box off
    if plotLoop < nPlots
        set(axPlots(plotLoop), 'XTickLabel', []);
    end
end
xlabel(axPlots(end), 'Time from block onset (s)')
linkaxes(axPlots, 'x');
xlim(axPlots(end), [blockTimeS(1) blockTimeS(end)]);

%% ---------------- Save ----------------
if saveBlocks
    [~, hekaBase] = fileparts(hekaFile);
    blocksFile = fullfile(pwd, [hekaBase '_laserBlocks.mat']);
    save(blocksFile, 'blockMats', 'blockTimeS', 'blockOnIdx', 'blockOffIdx', ...
         'blockDurS', 'normaliseTo405', 'preBlockSec', 'postBlockSec', ...
         'baselineWinSec', 'hekaFile');
    fprintf('Saved %s\n', blocksFile);
end
