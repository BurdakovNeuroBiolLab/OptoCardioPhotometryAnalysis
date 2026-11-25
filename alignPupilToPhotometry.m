%% alignPupilToPhotometry
% Puts the pupil diameter (one value per video frame, from pupilCircleFit)
% on the HEKA time base and adds it as column 11 of extractedMat, next to
% the photometry traces. Also builds extractedMatCrop, cropped to the
% camera period.
%
% Same logic as alignPupilToHeka.m in OptoCardioInscopixAnalysis:
% the camera sync (Cam Synch, column 10) runs from its first rising edge to
% its last falling edge for exactly the period covered by the pupil CSV.
% The video frames are spread evenly over that period; every HEKA sample
% takes the diameter of the video frame it falls in (NaN outside the
% camera period and in frames removed by pupilCircleFit, e.g. blinks).
%
% Needs in the workspace
%   from A_Extract_HekaData: extractedMat (nSamples x 10), colTitles,
%       sampleIntervalMs
%   pupilDiameter   cleaned diameter per video frame (pixels), nPupil x 1
%
% Adds
%   extractedMat(:, 11)  pupil diameter (px), colTitles{11}
%   extractedMatCrop     extractedMat cropped to camRise:camFall
%   camRise, camFall     first and last HEKA sample of the camera period
%   hekaMs               HEKA time of every sample (ms, 0 = first sample)
%   pupilHekaT           centre of each video frame in HEKA time (ms)
%   pupilPeriodMs        video frame period implied by the camera sync (ms)

camSynchCol = 10;
pupilCol    = 11;

if size(extractedMat, 2) >= pupilCol
    error('extractedMat already has a pupil column. Re-run A_Extract_HekaData first.');
end

%% Camera period from the Cam Synch channel
[camRise, camFall] = syncOnOff(extractedMat(:, camSynchCol), 'Cam Synch');

nHekaSamples = size(extractedMat, 1);
hekaMs   = (0:nHekaSamples-1)' * sampleIntervalMs;
camOnMs  = hekaMs(camRise);
camOffMs = hekaMs(camFall) + sampleIntervalMs;     % end of the last high sample

%% Spread the video frames evenly over the camera period
pupilDiameter = pupilDiameter(:);
nPupil        = numel(pupilDiameter);
pupilPeriodMs = (camOffMs - camOnMs) / nPupil;
pupilHekaT    = camOnMs + ((1:nPupil)' - 0.5) * pupilPeriodMs;

pupilFrameIdx = floor((hekaMs - camOnMs) / pupilPeriodMs) + 1;
inCam         = pupilFrameIdx >= 1 & pupilFrameIdx <= nPupil;
pupilOnHeka   = nan(nHekaSamples, 1);
pupilOnHeka(inCam) = pupilDiameter(pupilFrameIdx(inCam));

extractedMat(:, pupilCol) = pupilOnHeka;
colTitles{pupilCol}       = 'Pupil diameter';
extractedMatCrop          = extractedMat(camRise:camFall, :);

%% Diagnostic: if the camera sync pulses once per frame, the pulse count
% should match the number of rows in the pupil CSV.
camPulses = nnz(diff(extractedMat(camRise:camFall, camSynchCol)) == 1) + 1;
fprintf(['Pupil: %d video frames over %.2f s of camera sync -> %.2f Hz ' ...
         '(%d camera sync pulses)\n'], ...
        nPupil, (camOffMs - camOnMs) / 1000, 1000 / pupilPeriodMs, camPulses);

%% ---------------------------------------------------------------------
function [riseIdx, fallIdx] = syncOnOff(syncSig, sigName)
% First sample of the first rising edge and last high sample before the
% last falling edge of a binarised sync signal.
syncEdges = diff(syncSig);
riseIdx   = find(syncEdges == 1, 1, 'first') + 1;
fallIdx   = find(syncEdges == -1, 1, 'last');
if isempty(riseIdx) || isempty(fallIdx)
    error('%s has no complete rise and fall', sigName);
end
end
