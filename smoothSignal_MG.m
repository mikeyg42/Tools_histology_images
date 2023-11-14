function smoothSignal = smoothSignal_MG(counts)

signal = double(counts);

for m = 1:4
    vals = [];
    locsPeaks = [];
    locsValleys = [];
    pks = [];

    if length(signal) >= 3
        [pks, locsPeaks] = findpeaks(signal);
        [vals, locsValleys] = findpeaks(-signal);
    else
        error('Data set must contain at least 3 samples.');
    end

    vals = -vals; % Correct the sign for the valleys

    resampledSignalX = linspace(1, length(signal), length(signal));

    % Include the first and last points
    keyPointsXa = [1, locsPeaks, length(signal)];
    keyPointsYa = [signal(1), pks, signal(end)];
    if any(~isfinite(keyPointsXa)) || any(diff(keyPointsXa) <= 0)
        error('x must be finite and strictly increasing.');
    end
    resampledSignalYa = interp1(keyPointsXa, keyPointsYa, resampledSignalX, 'pchip');

    keyPointsXb = [1, locsValleys, length(signal)];
    keyPointsYb = [signal(1), vals, signal(end)];
    if any(~isfinite(keyPointsXb)) || any(diff(keyPointsXb) <= 0)
        error('x must be finite and strictly increasing.');
    end
    resampledSignalYb = interp1(keyPointsXb, keyPointsYb, resampledSignalX, 'pchip');

    resampledSignalY = (resampledSignalYa+resampledSignalYb)/2;

    if m ==1
        signal = resampledSignalY;
        continue
    end

    [a, ~] = findpeaks(resampledSignalY);

    if length(a) < 2
        smoothSignal = signal; % revert to previous iteration
        break

    elseif length(a) < 10
        smoothSignal = resampledSignalY;
        break

    else
        signal = resampledSignalY;

    end % end if block
end % ends the for loop
end