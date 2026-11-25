function alpha1 = dfa_alpha1_local(ibiSeries)
% DFA_ALPHA1_LOCAL: short-term DFA exponent for an IBI time series.
% Same as dfa_alpha1_local.m in OptoCardioInscopixAnalysis.
% Box sizes run from 4 beats up to floor(N/4), capped at 16 beats (the
% standard alpha1 range), so a long series still gives alpha1 rather than
% a mix of short- and long-term scaling. Needs >= 8 values.

    ibiSeries = ibiSeries(:);
    nValues   = numel(ibiSeries);

    % Center and integrate (ibiProfile)
    ibiProfile = cumsum(ibiSeries - mean(ibiSeries));

    % Box sizes: minimum 4, maximum floor(N/4), at most 16
    boxMin   = 4;
    boxMax   = min(floor(nValues/4), 16);
    boxSizes = boxMin:boxMax;
    fluct    = nan(size(boxSizes));

    for boxLoop = 1:numel(boxSizes)
        boxSize = boxSizes(boxLoop);
        nSeg    = floor(nValues / boxSize);
        if nSeg < 2
            continue;
        end

        rmsSeg  = zeros(nSeg, 1);
        boxTime = (1:boxSize)';
        for segLoop = 1:nSeg
            segIdx   = (segLoop-1)*boxSize + (1:boxSize);
            trendFit = polyfit(boxTime, ibiProfile(segIdx), 1);    % linear trend
            rmsSeg(segLoop) = sqrt(mean((ibiProfile(segIdx) - polyval(trendFit, boxTime)).^2));
        end
        fluct(boxLoop) = sqrt(mean(rmsSeg.^2));
    end

    validBox = ~isnan(fluct) & fluct > 0;
    if sum(validBox) < 2
        alpha1 = NaN;
        return;
    end

    slopeFit = polyfit(log(boxSizes(validBox))', log(fluct(validBox))', 1);
    alpha1   = slopeFit(1);
end
