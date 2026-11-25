
%% A_Extract_HekaData
% Loads HEKA-exported sweeps (*L.mat), concatenates them per channel,
% binarises the rotary, laser and camera-sync channels and demultiplexes
% the two interleaved 405/465 nm photometry channels.
%
% Output: extractedMat (nSamples x 10)
%   Col 1  : Heart rate (raw)
%   Col 2  : Rotary 1 (binarised)
%   Col 3  : Laser (binarised)
%   Col 4  : Phot 1 raw
%   Col 5  : Phot 1 405 (demultiplexed, spikes removed)
%   Col 6  : Phot 1 465 (demultiplexed, spikes removed)
%   Col 7  : Phot 2 raw
%   Col 8  : Phot 2 405 (demultiplexed, spikes removed)
%   Col 9  : Phot 2 465 (demultiplexed, spikes removed)
%   Col 10 : Cam sync (binarised)
% Last updated 25.11.2026

%% Settings
nChannels = 9;         % channels per sweep in the HEKA export
sampleIntervalMs = 2.5;  % HEKA sample interval (ms)

senseCheckGraphs = 1;  % 1 = plot raw vs demultiplexed photometry
cutTheFront      = 0;  % 1 = drop the first cutTheFrontVal samples
cutTheFrontVal   = 5000;

limitLength = 0;       % 1 = keep only the first superLength samples
superLength = 100000;

% Channel numbers in the HEKA export (columns of superRigMatMaster)
heartRateChan = 1;
rotary1Chan   = 2;
laserChan     = 3;
camSyncChan   = 8;

% Photometry demultiplexing: one LED cycle is 40 samples long. These are
% the sample positions within the cycle where each LED's signal is read.
samplesPerCycle = 40;
phot1_405start  = 30;  % CHECK THIS
phot1_465start  = 10;
phot2_405start  = 30;
phot2_465start  = 10;

spikeSDThresh = 10;    % samples above median + this many SDs become NaN

%% Find files, skipping macOS '._' resource files
fileNames = dir('*L.mat');
fileNames = fileNames(~startsWith({fileNames.name}, '._'));
numFiles  = numel(fileNames);

filesToProcess = 1;    % change to 1:numFiles to process every file

for fileLoop = filesToProcess

    fileName = fileNames(fileLoop).name;
    fprintf('Processing file %d: %s\n', fileLoop, fileName);

    %% Photometry channel numbers for this recording type
    if contains(fileName, '_t45L')
        photTrace1_Num = 4;
        photTrace2_Num = 5;
    else
        error('No photometry channel mapping defined for %s', fileName);
    end

    %% Load and drop the timestamp column of every sweep (2.5 ms/sample)
    rawCell = struct2cell(load(fileName));
    rawCell = cellfun(@(sweepData) sweepData(:, 2:end), rawCell, 'UniformOutput', false);

    % Remove the last sweep (one cell per channel) if it is incomplete
    if size(rawCell{end}, 1) ~= size(rawCell{1}, 1)
        rawCell(end-nChannels+1:end) = [];
    end

    %% Concatenate sweeps into one column per channel
    channelCols = cell(1, nChannels);
    for channelLoop = 1:nChannels
        channelCols{channelLoop} = vertcat(rawCell{channelLoop:nChannels:end});
    end
    superRigMatMaster = [channelCols{:}];

    if limitLength
        superRigMatMaster = superRigMatMaster(1:superLength, :);
    end

    if contains(fileName, 'absLevelCheck') || cutTheFront
        disp('Cutting first section')
        superRigMatMaster(1:cutTheFrontVal, :) = [];
    end

    nSamples = size(superRigMatMaster, 1);

    %% Figure 1: all raw HEKA traces
    figure('Renderer', 'painters', 'Position', [0 350 400 700]);
    for channelLoop = 1:nChannels
        subplot(nChannels, 1, channelLoop);
        plot(superRigMatMaster(:, channelLoop), '-k')
    end

    %% Photometry: demultiplex 405/465 and remove noise spikes
    phot1Trace = superRigMatMaster(:, photTrace1_Num);
    phot2Trace = superRigMatMaster(:, photTrace2_Num);

    phot1_405 = removeSpikes(demultiplexPhotometry(phot1Trace, phot1_405start, samplesPerCycle), spikeSDThresh);
    phot1_465 = removeSpikes(demultiplexPhotometry(phot1Trace, phot1_465start, samplesPerCycle), spikeSDThresh);
    phot2_405 = removeSpikes(demultiplexPhotometry(phot2Trace, phot2_405start, samplesPerCycle), spikeSDThresh);
    phot2_465 = removeSpikes(demultiplexPhotometry(phot2Trace, phot2_465start, samplesPerCycle), spikeSDThresh);

    %% Build the output matrix
    extractedMat = zeros(nSamples, 10);
    extractedMat(:, 1)  = superRigMatMaster(:, heartRateChan);
    extractedMat(:, 2)  = binariseTrace(superRigMatMaster(:, rotary1Chan));
    extractedMat(:, 3)  = binariseTrace(superRigMatMaster(:, laserChan));
    extractedMat(:, 4)  = phot1Trace;
    extractedMat(:, 5)  = phot1_405;
    extractedMat(:, 6)  = phot1_465;
    extractedMat(:, 7)  = phot2Trace;
    extractedMat(:, 8)  = phot2_405;
    extractedMat(:, 9)  = phot2_465;
    extractedMat(:, 10) = binariseTrace(superRigMatMaster(:, camSyncChan));

    colTitles = {'HR', 'Rotary 1', 'Laser', ...
                 'Phot1 raw', 'Phot1 405', 'Phot1 465', ...
                 'Phot2 raw', 'Phot2 405', 'Phot2 465', 'Cam sync'};

    %% Sense-check plots
    if senseCheckGraphs
        plotRange = 1:500;

        figure('Renderer', 'painters', 'Position', [0 0 400 300]);

        subplot(2,2,1)
        plot(phot1Trace(plotRange), '-b'); hold on;
        plot(phot1_405(plotRange), '-m')
        plot(phot1_465(plotRange), '-g')
        legend('raw','405','465')
        ylabel('phot 1 trace')

        subplot(2,2,2)
        plot(phot1Trace, '-ob'); hold on;
        plot(phot1_405, '-om')
        plot(phot1_465, '-og')

        subplot(2,2,3)
        plot(phot2Trace(plotRange), '-b'); hold on;
        plot(phot2_405(plotRange), '-m')
        plot(phot2_465(plotRange), '-g')
        legend('raw','405','465')
        ylabel('phot 2 trace')

        subplot(2,2,4)
        plot(phot2Trace, '-b'); hold on;
        plot(phot2_405, '-m')
        plot(phot2_465, '-g')
    end

    %% Overview plot of every column of extractedMat
    nCols = size(extractedMat, 2);
    figure('Renderer', 'painters', 'Position', [500 200 400 800]);
    for plotLoop = 1:nCols
        subplot(nCols, 1, plotLoop)
        plot(extractedMat(:, plotLoop))
        ylabel(colTitles{plotLoop})
    end

end

%% Local functions (must stay at the end of the script)

function binaryTrace = binariseTrace(trace)
% Returns 1 where the trace is above the midpoint between its min and max.
threshold   = min(trace) + (max(trace) - min(trace)) / 2;
binaryTrace = double(trace > threshold);
end

function demux = demultiplexPhotometry(trace, startIdx, period)
% Rebuilds one LED's signal from an interleaved photometry trace. Every
% `period` samples, the mean of the 3 samples around that LED's position
% is held for the whole cycle.
traceLength = numel(trace);
demux = zeros(traceLength, 1);
demux(1:startIdx+period) = trace(startIdx);
lastIdx = startIdx + period;
for cycleStart = startIdx+period : period : traceLength-period
    demux(cycleStart:cycleStart+period-1) = mean(trace(cycleStart-1:cycleStart+1));
    lastIdx = cycleStart + period - 1;
end
demux(lastIdx+1:end) = demux(lastIdx);   % hold last value to the end
end

function trace = removeSpikes(trace, nSD)
% Sets samples more than nSD standard deviations above the median to NaN.
trace(trace > median(trace) + nSD * std(trace)) = NaN;
end
