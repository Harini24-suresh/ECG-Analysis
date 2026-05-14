function benchmark_rpeak_detection(record_id, fs)
    % Load ECG
    dat_file = [record_id, '.dat'];
    fid = fopen(dat_file, 'rb');
    if fid == -1
        fprintf("❌ Could not open %s.dat\n", record_id);
        return;
    end
    ecg = fread(fid, inf, 'uint8');
    fclose(fid);
    ecg = double(ecg);
    ecg = ecg - mean(ecg);  % Remove DC

    % Run R-peak detection
    r_peaks = detect_rpeaks_simple(ecg, fs);

    % Load true annotations
    qrs_file = [record_id, '.qrs'];
    if ~isfile(qrs_file)
        fprintf("⚠️ No .qrs annotation available for record %s. Skipping benchmarking.\n", record_id);
        return;
    end
    true_peaks = load(qrs_file);  % must be plain text list of sample indices

    % Benchmark
    tolerance = round(0.1 * fs);  % 100ms window
    TP = 0; FP = 0; FN = 0;
    matched = false(size(true_peaks));

    for i = 1:length(r_peaks)
        diffs = abs(true_peaks - r_peaks(i));
        [min_diff, idx] = min(diffs);
        if min_diff <= tolerance && ~matched(idx)
            TP = TP + 1;
            matched(idx) = true;
        else
            FP = FP + 1;
        end
    end

    FN = sum(~matched);
    Sensitivity = TP / (TP + FN) * 100;
    Precision = TP / (TP + FP) * 100;
    F1 = 2 * (Sensitivity * Precision) / (Sensitivity + Precision);

    % Print Results
    fprintf("=== Benchmark for %s ===\n", record_id);
    fprintf("TP = %d | FP = %d | FN = %d\n", TP, FP, FN);
    fprintf("Sensitivity = %.2f%% | Precision = %.2f%% | F1 = %.2f%%\n", Sensitivity, Precision, F1);
end
