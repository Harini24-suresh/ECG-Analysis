function r_peaks = detect_rpeaks_simple(ecg, fs)
    window_size = round(0.150 * fs);
    ecg_filt = filter(ones(1, window_size)/window_size, 1, ecg);

    diff_ecg = diff(ecg_filt);
    squared = diff_ecg .^ 2;
    integrated = conv(squared, ones(1, window_size)/window_size, 'same');
    threshold = 0.6 * max(integrated);
    qrs_locs = find(integrated > threshold);

    min_dist = round(0.25 * fs);  % 250 ms refractory
    r_peaks = [];
    if ~isempty(qrs_locs)
        r_peaks(1) = qrs_locs(1);
        for i = 2:length(qrs_locs)
            if (qrs_locs(i) - r_peaks(end)) > min_dist
                r_peaks(end+1) = qrs_locs(i);
            end
        end
    end
end
