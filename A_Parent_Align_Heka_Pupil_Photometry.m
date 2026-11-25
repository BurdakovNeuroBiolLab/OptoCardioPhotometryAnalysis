%% A_Parent_Align_Heka_Pupil_Photometry
% Parent script: loads the pupil video data and the HEKA recording from the
% MATLAB current folder and puts the pupil diameter on the same time base
% as the photometry traces.
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
%   4. Figure of photometry, laser and pupil over the camera period
%   5. Save
%
% The current folder must contain exactly one of each:
%   *zoompup*.csv (DeepLabCut pupil) and the HEKA export *L.mat.
%   Files starting with ._ are ignored.
%
% Output (in the workspace, and saved to <HEKA name>_pupilAligned.mat)
%   extractedMat      (nSamples x 11), 2.5 ms per row
%   extractedMatCrop  extractedMat cropped to the camera period
%   Columns: 1 HR          5 Phot 1 405   9  Phot 2 465
%            2 Rotary 1    6 Phot 1 465   10 Cam synch
%            3 Laser       7 Phot 2 raw   11 Pupil diameter (px)
%            4 Phot 1 raw  8 Phot 2 405
%   colTitles         column names
%   hekaMs            time of every row of extractedMat (ms)
%   hekaMsCrop        time of every row of extractedMatCrop (ms)
%   pupilHekaT        centre of each video frame in HEKA time (ms)

clear; clc;

%% ---------------- Settings ----------------
expectedFps = [];     % camera frame rate (Hz). If set, warn when the rate
                      % implied by the Cam Synch period differs by more
                      % than fpsTol. [] = no check.
fpsTol      = 0.02;   % fractional tolerance for the frame-rate check
saveAligned = true;   % save the aligned data to a .mat in the current folder

%% ---------------- Put the child scripts on the path ----------------
alignRepoDir = fileparts(mfilename('fullpath'));
if ~isempty(alignRepoDir)
    addpath(alignRepoDir);
end
childScripts = {'pupilCircleFit', 'A_Extract_HekaData', 'alignPupilToPhotometry'};
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

if ~isempty(expectedFps)
    impliedFps = 1000 / pupilPeriodMs;
    if abs(impliedFps / expectedFps - 1) > fpsTol
        warning(['Cam Synch period implies %.2f Hz but expectedFps is %.2f Hz. ' ...
                 'Check that the pupil CSV and the HEKA file are from the same ' ...
                 'recording.'], impliedFps, expectedFps);
    end
end

%% ---------------- 4. Photometry, laser and pupil over the camera period
plotCols   = [6 5 9 8 3 11];
plotColour = {[0 0.6 0], [0.7 0 0.7], [0 0.6 0], [0.7 0 0.7], [0 0 0.8], [0.9 0.5 0]};
cropTimeS  = (hekaMsCrop - hekaMsCrop(1)) / 1000;    % s from camera start

nPlots = numel(plotCols);
figure('Name', 'Photometry and pupil', 'Renderer', 'painters', ...
       'Position', [900 100 700 800]);
axPlots = gobjects(nPlots, 1);
for plotLoop = 1:nPlots
    axPlots(plotLoop) = subplot(nPlots, 1, plotLoop);
    plot(cropTimeS, extractedMatCrop(:, plotCols(plotLoop)), ...
         'Color', plotColour{plotLoop})
    ylabel(colTitles{plotCols(plotLoop)})
    box off
    if plotLoop < nPlots
        set(axPlots(plotLoop), 'XTickLabel', []);
    end
end
xlabel(axPlots(end), 'Time from camera start (s)')
linkaxes(axPlots, 'x');
xlim(axPlots(end), [cropTimeS(1) cropTimeS(end)]);

%% ---------------- 5. Save ----------------
if saveAligned
    [~, hekaBase] = fileparts(hekaFile);
    alignedFile = fullfile(pwd, [hekaBase '_pupilAligned.mat']);
    save(alignedFile, 'extractedMat', 'extractedMatCrop', 'colTitles', ...
         'hekaMs', 'hekaMsCrop', 'camRise', 'camFall', 'pupilDiameter', ...
         'pupilHekaT', 'pupilPeriodMs', 'sampleIntervalMs', 'hekaFile', 'pupilFile');
    fprintf('Saved %s\n', alignedFile);
end
