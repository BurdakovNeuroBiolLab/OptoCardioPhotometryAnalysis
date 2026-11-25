function metrics = continuous_hrv_mouse(beatTimes, windowLengthSec, stepSec, ibiRangeSec, ibiRelTol)
% CONTINUOUS_HRV_MOUSE
%   Compute continuous HRV measures from mouse beat times in sliding windows.
%   Same as continuous_hrv_mouse.m in OptoCardioInscopixAnalysis.
%
% INPUTS
%   beatTimes       - vector of heartbeat (R-peak) times in seconds
%   windowLengthSec - length of sliding window in seconds
%   stepSec         - step between window starts in seconds
%   ibiRangeSec     - (optional) [min max] plausible inter-beat interval (s).
%                     IBIs outside it (missed / extra beats) are set to NaN
%                     and left out of every metric. Default: keep all.
%                     Mouse example: 60 ./ [900 250] = [0.067 0.24].
%   ibiRelTol       - (optional) also reject IBIs that differ from the running
%                     median of 21 IBIs by more than this fraction (e.g. 0.3).
%                     Catches missed beats (IBI doubles), which can still be
%                     inside ibiRangeSec. Default: off. See flag_bad_ibi.m.
%
% OUTPUT struct 'metrics' (all column vectors, one row per window)
%   metrics.t          - centre time of each window (s)
%   metrics.HR         - heart rate per window (bpm), 60 / mean IBI
%   metrics.RMSSD      - root mean square of successive IBI differences (s)
%   metrics.CVSD       - RMSSD / mean IBI
%   metrics.CV_IBI     - SDNN / mean IBI
%   metrics.alpha1_DFA - DFA alpha1 (unreliable for windows < ~60 beats)
%   metrics.logRMSSD   - log(RMSSD)
%   metrics.HRV_resid  - residual log(RMSSD) after regressing out HR
%   metrics.nIBI       - number of IBIs used per window
%
% Successive differences are only taken between two valid, consecutive IBIs,
% so a rejected IBI never contributes a difference.

    if nargin < 4 || isempty(ibiRangeSec)
        ibiRangeSec = [0 Inf];
    end
    if nargin < 5
        ibiRelTol = [];
    end
    beatTimes = beatTimes(:);

    % Flag bad IBIs once over the whole recording; IBI k is beat k -> k+1
    ibiAll = diff(beatTimes);
    badAll = flag_bad_ibi(ibiAll, ibiRangeSec, ibiRelTol, 21);

    % Window start times
    winStarts = beatTimes(1):stepSec:(beatTimes(end) - windowLengthSec);
    nWin = numel(winStarts);

    % Preallocate
    tCenter   = nan(nWin, 1);
    HR        = nan(nWin, 1);
    RMSSD     = nan(nWin, 1);
    CVSD      = nan(nWin, 1);
    CV_IBI    = nan(nWin, 1);
    alpha1    = nan(nWin, 1);
    logRMSSD  = nan(nWin, 1);
    nIBI      = zeros(nWin, 1);

    for winLoop = 1:nWin
        winStart = winStarts(winLoop);
        winEnd   = winStart + windowLengthSec;
        tCenter(winLoop) = winStart + windowLengthSec / 2;

        % Beats within the moving window
        winBeats = find(beatTimes > winStart & beatTimes <= winEnd);
        if numel(winBeats) < 2
            continue
        end
        IBI = ibiAll(winBeats(1):winBeats(end) - 1);          % inter-beat intervals (s)
        IBI(badAll(winBeats(1):winBeats(end) - 1)) = NaN;
        goodIBI = IBI(~isnan(IBI));
        nIBI(winLoop) = numel(goodIBI);
        if isempty(goodIBI)
            continue
        end

        meanIBI = mean(goodIBI);
        HR(winLoop) = 60 / meanIBI;             % bpm

        diffIBI = diff(IBI);                    % NaN next to a rejected IBI
        diffIBI = diffIBI(~isnan(diffIBI));
        if ~isempty(diffIBI)
            RMSSD(winLoop) = sqrt(mean(diffIBI .^ 2));
        end
        sdnn = sqrt(mean((goodIBI - meanIBI) .^ 2));   % population SD, as std(x, 1)

        CVSD(winLoop)     = RMSSD(winLoop) / meanIBI;
        CV_IBI(winLoop)   = sdnn / meanIBI;
        alpha1(winLoop)   = dfa_alpha1_local(goodIBI);
        logRMSSD(winLoop) = log(RMSSD(winLoop));
    end

    % Residual HRV: regress logRMSSD on HR across windows
    validWin  = isfinite(HR) & isfinite(logRMSSD);     % excludes log(0) = -Inf
    HRV_resid = nan(size(HR));
    if sum(validWin) >= 3
        designMat = [ones(sum(validWin), 1), HR(validWin)];
        regCoeffs = designMat \ logRMSSD(validWin);
        HRV_resid(validWin) = logRMSSD(validWin) - designMat * regCoeffs;
    end

    % Output
    metrics = struct();
    metrics.t          = tCenter;
    metrics.HR         = HR;
    metrics.CVSD       = CVSD;
    metrics.CV_IBI     = CV_IBI;
    metrics.RMSSD      = RMSSD;
    metrics.alpha1_DFA = alpha1;
    metrics.logRMSSD   = logRMSSD;
    metrics.HRV_resid  = HRV_resid;
    metrics.nIBI       = nIBI;
end
