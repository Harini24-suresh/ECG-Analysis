function benchmark_rpeak_detection(record_id, fs)

    fprintf("\n=== RECORD: %s ===\n", record_id);

    % Load ECG
    fid = fopen([record_id, '.dat'], 'rb');
    if fid == -1
        error('Cannot open ECG file: %s.dat', record_id);
    end
    ecg = double(fread(fid, inf, 'uint8'));
    fclose(fid);
    ecg = ecg - mean(ecg);

    % Run your Pan-Tompkins algorithm (replace with your function)
    detected_peaks = pan_tompkins(ecg, fs);

    % Load true QRS annotations
    try
        true_peaks = load([record_id, '.qrs']);  % assumed to be sample indices
    catch
        fprintf("⚠️ No .qrs file for %s. Skipping benchmarking.\n", record_id);
        return;
    end

    % Match with tolerance
    tolerance = round(0.05 * fs);  % 50 ms
    TP = 0; FP = 0; FN = 0;

    matched_detected = false(size(detected_peaks));
    matched_true = false(size(true_peaks));

    % True Positives and False Negatives
    for i = 1:length(true_peaks)
        diffs = abs(detected_peaks - true_peaks(i));
        [min_diff, idx] = min(diffs);
        if min_diff <= tolerance
            TP = TP + 1;
            matched_detected(idx) = true;
            matched_true(i) = true;
        else
            FN = FN + 1;
        end
    end

    % False Positives (unmatched detections)
    FP = sum(~matched_detected);

    % Metrics
    precision = TP / (TP + FP);
    recall = TP / (TP + FN);
    F1 = 2 * (precision * recall) / (precision + recall);

    % Display results
    fprintf("✅ TP = %d | ❌ FP = %d | ❗ FN = %d\n", TP, FP, FN);
    fprintf("📊 Precision = %.2f%% | Recall = %.2f%% | F1 = %.2f%%\n", ...
            precision*100, recall*100, F1*100);
end
