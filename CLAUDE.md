# CLAUDE.md

Instructions for Claude when working in this repository.

## Code style

### Meaningful variable names — no single letters

Do not use single-letter or meaningless names such as `i`, `ii`, `j`, `k`, `x`, `n`
for loop counters, function arguments, anonymous-function inputs or any other
variable. They are hard for a reader to decipher at a glance. Use a name that
says what the variable is.

Loop counters are named after what they loop over, ending in `Loop`:

```matlab
% Not this
for k = 1:numel(plotCols)
    subplot(numel(plotCols), 1, k)
    plot(binarisedMat(:, plotCols(k)))
    xlabel(xLabTitles{k})
end

% This
for plotLoop = 1:numel(plotCols)
    subplot(numel(plotCols), 1, plotLoop)
    plot(binarisedMat(:, plotCols(plotLoop)))
    xlabel(xLabTitles{plotLoop})
end
```

Other examples: `channelLoop`, `fileLoop`, `trialLoop`; `trace` instead of `x`;
`traceLength` instead of `L`; `@(sweepData) ...` instead of `@(x) ...`.
