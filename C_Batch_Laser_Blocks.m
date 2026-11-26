%% C_Batch_Laser_Blocks
% Runs A_Parent_Align_Heka_Pupil_Photometry and B_Align_Laser_Blocks in
% every day folder of every mouse, then averages the laser-block traces
% within each mouse and across mice, separately for ctrl and opto days.
%
% Run from the main folder. Expected layout:
%   main folder / m947 / d1 opto
%                      / d2 ctrl
%               / m946 / d1 ctrl
%                      / d2 opto
% The treatment is taken from the day folder name (it must contain 'opto'
% or 'ctrl'); other folders are skipped. Each day folder needs one HEKA
% export (*L.mat) and one pupil CSV (*zoompup*.csv).
%
% Settings for the HEKA / heart rate / pupil processing are in
% A_Parent_Align_Heka_Pupil_Photometry, and for the laser blocks
% (normaliseTo405, zscorePupil, windows, baseline) in B_Align_Laser_Blocks.
% They apply to every day.
%
% Output in the main folder, for each treatment (ctrl / opto):
%   laserBlocks_<treatment>.mat
%     groupMats       struct, one nMice x nTime matrix per signal (phot1_465z,
%                     phot2_465z, HR, RMSSD, pupil, run, laser): every row
%                     is one mouse, averaged over all blocks of all its
%                     days of that treatment
%     groupTimeS      1 x nTime, time from block onset (s)
%     mouseNames      mouse of every row of groupMats
%     mouseDays       day folders averaged into every row
%     pooledBlocks    struct, every block of every day, one row per block
%     pooledMouseIdx  row of groupMats (mouse) of every pooled block
%     plotNames       y-axis labels
%   laserBlocks_<treatment>.fig / .png
%     same layout as the B_Align_Laser_Blocks figure: one colour per
%     mouse, mean across mice thick black, median thick dotted black

clear; clc;

%% ---------------- Settings ----------------
mainFolder      = pwd;
treatments      = {'ctrl', 'opto'};
closeDayFigures = true;   % close each day's figures once they are saved
saveGroupFigs   = true;   % save the ctrl / opto figures as .fig and .png

%% ---------------- Put the repository on the path ----------------
repoDir = fileparts(mfilename('fullpath'));
if ~isempty(repoDir)
    addpath(repoDir);
end

%% ---------------- Run every day folder ----------------
mouseDirs = dir(mainFolder);
mouseDirs = mouseDirs([mouseDirs.isdir] & ~startsWith({mouseDirs.name}, '.'));

dayResults = struct('mouseName', {}, 'dayName', {}, 'treatment', {}, ...
                    'blockMats', {}, 'blockTimeS', {}, 'blockDurS', {}, ...
                    'blockPlotFields', {}, 'blockPlotNames', {}, ...
                    'blocksFromTemplate', {});
failedDays = {};

for mouseLoop = 1:numel(mouseDirs)
    mouseName = mouseDirs(mouseLoop).name;
    dayDirs   = dir(fullfile(mainFolder, mouseName));
    dayDirs   = dayDirs([dayDirs.isdir] & ~startsWith({dayDirs.name}, '.'));

    for dayLoop = 1:numel(dayDirs)
        dayName = dayDirs(dayLoop).name;
        if contains(lower(dayName), 'opto')
            treatment = 'opto';
        elseif contains(lower(dayName), 'ctrl')
            treatment = 'ctrl';
        else
            fprintf('Skipping %s / %s (no ''opto'' or ''ctrl'' in the name)\n', mouseName, dayName);
            continue
        end

        fprintf('\n========== %s %s ==========\n', mouseName, dayName);
        try
            dayOut = runOneDay(fullfile(mainFolder, mouseName, dayName));
            dayOut.mouseName = mouseName;
            dayOut.dayName   = dayName;
            dayOut.treatment = treatment;
            dayResults(end + 1) = orderfields(dayOut, dayResults); %#ok<SAGROW>
        catch dayError
            warning('%s %s failed: %s', mouseName, dayName, dayError.message);
            failedDays{end + 1} = sprintf('%s %s: %s', mouseName, dayName, dayError.message); %#ok<SAGROW>
        end
        cd(mainFolder);
        if closeDayFigures
            close all
        end
    end
end

%% ---------------- Average within mouse, then across mice ----------------
for treatmentLoop = 1:numel(treatments)
    treatment = treatments{treatmentLoop};
    isTreat   = find(strcmp({dayResults.treatment}, treatment));
    if isempty(isTreat)
        warning('No %s days were processed.', treatment);
        continue
    end

    plotFields = dayResults(isTreat(1)).blockPlotFields;
    plotNames  = dayResults(isTreat(1)).blockPlotNames;
    allFields  = fieldnames(dayResults(isTreat(1)).blockMats);

    % Common time axis: every day's window starts at -preBlockSec with the
    % same sample interval; shorter windows are padded with NaN.
    timeLengths = arrayfun(@(dayIdx) numel(dayResults(dayIdx).blockTimeS), isTreat);
    [nTime, longestDay] = max(timeLengths);
    groupTimeS = dayResults(isTreat(longestDay)).blockTimeS;
    laserDurS  = max(arrayfun(@(dayIdx) max(dayResults(dayIdx).blockDurS), isTreat));

    mouseNames = unique({dayResults(isTreat).mouseName}, 'stable');
    nMice      = numel(mouseNames);
    mouseDays  = cell(nMice, 1);

    groupMats      = struct();
    pooledBlocks   = struct();
    pooledMouseIdx = [];
    for fieldLoop = 1:numel(allFields)
        groupMats.(allFields{fieldLoop})    = nan(nMice, nTime);
        pooledBlocks.(allFields{fieldLoop}) = nan(0, nTime);
    end

    for mouseLoop = 1:nMice
        mouseDayIdx = isTreat(strcmp({dayResults(isTreat).mouseName}, mouseNames{mouseLoop}));
        mouseDays{mouseLoop} = {dayResults(mouseDayIdx).dayName};

        for fieldLoop = 1:numel(allFields)
            fieldName   = allFields{fieldLoop};
            mouseBlocks = nan(0, nTime);
            for dayIdx = mouseDayIdx
                dayMat = dayResults(dayIdx).blockMats.(fieldName);
                padded = nan(size(dayMat, 1), nTime);
                padded(:, 1:size(dayMat, 2)) = dayMat;
                mouseBlocks = [mouseBlocks; padded]; %#ok<AGROW>
            end
            groupMats.(fieldName)(mouseLoop, :) = mean(mouseBlocks, 1, 'omitnan');
            pooledBlocks.(fieldName) = [pooledBlocks.(fieldName); mouseBlocks];
        end
        nMouseBlocks   = size(mouseBlocks, 1);
        pooledMouseIdx = [pooledMouseIdx; repmat(mouseLoop, nMouseBlocks, 1)]; %#ok<AGROW>
    end

    groupTitle = sprintf('%s (n = %d mice: %s)', treatment, nMice, strjoin(mouseNames, ', '));
    groupFig   = plotLaserBlockFigure(groupTimeS, groupMats, plotFields, plotNames, ...
                                      laserDurS, groupTitle, true);

    groupBase = fullfile(mainFolder, ['laserBlocks_' treatment]);
    save([groupBase '.mat'], 'groupMats', 'groupTimeS', 'mouseNames', 'mouseDays', ...
         'pooledBlocks', 'pooledMouseIdx', 'plotFields', 'plotNames', 'laserDurS');
    fprintf('Saved %s.mat\n', groupBase);
    if saveGroupFigs
        savefig(groupFig, [groupBase '.fig']);
        print(groupFig, [groupBase '.png'], '-dpng', '-r300');
        fprintf('Saved %s.fig and .png\n', groupBase);
    end
end

%% ---------------- Summary ----------------
fprintf('\nProcessed %d day folders.\n', numel(dayResults));
templateDays = dayResults([dayResults.blocksFromTemplate]);
for templateLoop = 1:numel(templateDays)
    fprintf('  No laser found, template block times used: %s %s\n', ...
            templateDays(templateLoop).mouseName, templateDays(templateLoop).dayName);
end
if ~isempty(failedDays)
    fprintf('%d day folders failed:\n', numel(failedDays));
    fprintf('  %s\n', failedDays{:});
end

%% ---------------------------------------------------------------------
function dayOut = runOneDay(dayPath)
% Runs A and B in one day folder. Running them inside this function gives
% every day a fresh, empty workspace.
cd(dayPath);
batchMode = true; %#ok<NASGU>  read by A_Parent_Align_Heka_Pupil_Photometry
A_Parent_Align_Heka_Pupil_Photometry;
B_Align_Laser_Blocks;

dayOut = struct();
dayOut.blockMats          = blockMats;
dayOut.blockTimeS         = blockTimeS;
dayOut.blockDurS          = blockDurS;
dayOut.blockPlotFields    = blockPlotFields;
dayOut.blockPlotNames     = blockPlotNames;
dayOut.blocksFromTemplate = blocksFromTemplate;
end
