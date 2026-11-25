function bad = flag_bad_ibi(ibi, rangeSec, relTol, medWin)
% FLAG_BAD_IBI  Mark inter-beat intervals that come from missed or extra beats.
% Same as flag_bad_ibi.m in OptoCardioInscopixAnalysis.
%
%   bad = flag_bad_ibi(ibi, rangeSec, relTol, medWin)
%
%   ibi      - inter-beat intervals (s)
%   rangeSec - [min max] plausible IBI (s), e.g. 60 ./ [900 250]
%   relTol   - reject IBIs differing from the local median by more than this
%              fraction (e.g. 0.3). A missed beat doubles the IBI, an extra
%              beat roughly halves it, and both can still lie inside rangeSec.
%   medWin   - number of IBIs in the running median (e.g. 21)
%
%   bad      - logical, true for IBIs to leave out

    ibi = ibi(:);
    bad = ibi < rangeSec(1) | ibi > rangeSec(2) | isnan(ibi);
    if nargin >= 3 && ~isempty(relTol) && numel(ibi) >= 3
        ibiGood = ibi;
        ibiGood(bad) = NaN;
        ibiMedian = movmedian(ibiGood, medWin, 'omitnan');
        bad = bad | abs(ibi ./ ibiMedian - 1) > relTol;
    end
end
