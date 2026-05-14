% Cross-validation for R-peak detection across segments
clear; clc;

fs = 100;  % Sampling freq
records = {'a01', 'a02', 'a04', 'a19', 'b04', 'c04'};

for i = 1:length(records)
    record_id = records{i};
    fprintf("\n========= RECORD: %s =========\n", upper(record_id));

    % Load ECG
    try
        fid = fopen([record_id, '.dat'], 'rb');
        ecg = fread(fid, inf, 'uint8'); fclose(fid);
        ecg = double(ecg);
    catch
        fprintf("❌ Could not load ECG for %s\n", record_id);
        continue;
    end

    % Normalize
    ecg = ecg - mean(ecg);
    ecg = ecg / std(ecg);

    % Split into 3 segments
    L = length(ecg);
    segments = {
        ecg(1 : floor(L/3)), ...
        ecg(floor(L/3)+1 : floor(2*L/3)), ...
        ecg(floor(2*L/3)+1 : end)
    };

    % Apply your Pan-Tompkins QRS detection to each segment
    rpeaks = cell(1,3);
    counts = zeros(1,3);

    for j = 1:3
        seg = segments{j};

        % --- Filtering and detection (same as your main QRS logic) ---
        window_size = round(0.150 * fs);
        filt = filter(ones(1, window_size)/window_size, 1, seg);

        diff_ecg = diff(filt);
        squared = diff_ecg .^ 2;
        integrated = conv(squared, ones(1, window_size)/window_size, 'same');
        threshold = 0.6 * max(integrated);
        locs = find(integrated > threshold);

        min_dist = round(0.25 * fs);
        qrs = [];
        if ~isempty(locs)
            qrs(1) = locs(1);
            for k = 2:length(locs)
                if (locs(k) - qrs(end)) > min_dist
                    qrs(end+1) = locs(k);
                end
            end
        end

        rpeaks{j} = qrs;
        counts(j) = length(qrs);
    end

    % Print results
    fprintf("Segment 1: R-peaks = %d\n", counts(1));
    fprintf("Segment 2: R-peaks = %d\n", counts(2));
    fprintf("Segment 3: R-peaks = %d\n", counts(3));

    % Consistency Check
    if std(counts) > 0.2 * mean(counts)
        fprintf("❌ Inconsistent R-peak detection — Check signal quality\n");
    else
        fprintf("✅ Consistent R-peak detection\n");
    end
end

