% Multi-ECG Analysis Script (Apnea, HRV, Arrhythmia Detection with Subplots)
clear; clc; close all;

% Setup
fs = 100;  % Sampling frequency in Hz
records = {'a01', 'a02', 'a04','a19','b04','c04'};

% Results table
results = table('Size', [0 5], ...
    'VariableTypes', {'string','double','double','double','double'}, ...
    'VariableNames', {'Record','MeanHR','SDNN','RMSSD','IrregularCount'});

for r = 1:length(records)
    record_id = records{r};
    fprintf('\n=== Processing %s ===\n', record_id);

    %% Load ECG
    dat_file = [record_id, '.dat'];
    fid = fopen(dat_file, 'rb');
    if fid == -1
        warning('Cannot open file: %s', dat_file);
        continue;
    end
    ecg = fread(fid, inf, 'uint8');
    fclose(fid);
    ecg = double(ecg);
    ecg = ecg - mean(ecg);  % DC offset removal

    %% Load Apnea Annotations
    apn_file = [record_id, '.apn'];
    fid = fopen(apn_file, 'r');
    if fid == -1
        warning('Cannot open file: %s', apn_file);
        continue;
    end
    apnea_ann = fscanf(fid, '%d');
    fclose(fid);

    samples_per_min = fs * 60;
    apnea_mask = zeros(size(ecg));
    for i = 1:length(apnea_ann)
        if apnea_ann(i) == 1
            idx_start = (i - 1) * samples_per_min + 1;
            idx_end = min(i * samples_per_min, length(ecg));
            apnea_mask(idx_start:idx_end) = 1;
        end
    end

    %% Filtering
    window_size = round(0.150 * fs);
    ecg_filt = filter(ones(1, window_size)/window_size, 1, ecg);

    if strcmp(record_id, 'a19')
    t = (0:length(ecg)-1)/fs;
    figure('Name','ECG Signal Quality - a19','NumberTitle','off');
    plot(t, ecg, 'b');
    title('Raw ECG Signal - a19 (Check Signal Quality)');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

    %% QRS Detection
    diff_ecg = diff(ecg_filt);
    squared = diff_ecg .^ 2;
    integrated = conv(squared, ones(1, window_size)/window_size, 'same');
    threshold = 0.6 * max(integrated);
    qrs_locs = find(integrated > threshold);

    min_dist = round(0.25 * fs);  % 250ms refractory
    qrs_peaks = [];
    
    if ~isempty(qrs_locs)
        qrs_peaks(1) = qrs_locs(1);
        for i = 2:length(qrs_locs)
            if (qrs_locs(i) - qrs_peaks(end)) > min_dist
                qrs_peaks(end+1) = qrs_locs(i); %#ok<AGROW>
            end
        end
    end


    % === SAVE QRS PEAKS TO FILE ===
qrs_filename = [record_id, '.qrs'];
save(qrs_filename, 'qrs_peaks', '-ascii');
fprintf('[INFO] Saved R-peaks to %s\n', qrs_filename);

    %% HRV Metrics
    rr_intervals = diff(qrs_peaks) / fs;
    mean_rr = mean(rr_intervals);
    sdnn = std(rr_intervals);
    rmssd = sqrt(mean(diff(rr_intervals).^2));
    heart_rates = 60 ./ rr_intervals;
    rr_times = qrs_peaks(2:end) / fs;

    %% Arrhythmia Detection
    brady_idx = find(heart_rates < 60);
    tachy_idx = find(heart_rates > 100);
    irregular_idx = find(abs(diff(rr_intervals)) > 0.2);

    %% Store Results
    results = [results; {record_id, mean(heart_rates), sdnn, rmssd, length(irregular_idx)}];

    %% PLOTTING (All-in-One Figure for Each Record)
    t = (0:length(ecg)-1) / fs;
    figure('Name', ['Multi-Plot for Record ', record_id], 'NumberTitle', 'off');
    
    % ---- Subplot 1: ECG + Apnea ----
    subplot(3,1,1);
    plot(t, ecg_filt, 'b'); hold on;
    plot(t, max(ecg_filt) * apnea_mask, 'r', 'LineWidth', 1.2);
    title(['ECG with Apnea Overlay - ', record_id]);
    xlabel('Time (s)'); ylabel('Amplitude');
    legend('ECG','Apnea'); grid on;

    % ---- Subplot 2: ECG + QRS ----
    subplot(3,1,2);
    plot(t, ecg_filt); hold on;
    plot(qrs_peaks/fs, ecg_filt(qrs_peaks), 'ro');
    title(['QRS Detection - ', record_id]);
    xlabel('Time (s)'); ylabel('Amplitude');
    legend('ECG','QRS Peaks'); grid on;

    % ---- Subplot 3: Heart Rate ----
    subplot(3,1,3);
    plot(rr_times, heart_rates, 'b'); hold on;
    plot(rr_times(brady_idx), heart_rates(brady_idx), 'go');
    plot(rr_times(tachy_idx), heart_rates(tachy_idx), 'ro');
    plot(rr_times(irregular_idx+1), heart_rates(irregular_idx+1), 'mo');
    yline(mean(heart_rates), '--k', 'Mean HR');
    title(['Heart Rate with Arrhythmias - ', record_id]);
    xlabel('Time (s)'); ylabel('BPM');
    legend('HR','Brady (<60)','Tachy (>100)','Irregular','Mean HR');
    grid on;
end

%% Final Results
disp('==== HRV Results Table ====');
disp(results);

%% Compute 95% Confidence Intervals for SDNN and RMSSD

% Extract metrics
sdnn_vals = results.SDNN;
rmssd_vals = results.RMSSD;

% Sample size
n = length(sdnn_vals);
t_score = tinv(0.975, n-1);  % For 95% confidence

% SDNN CI
sdnn_mean = mean(sdnn_vals);
sdnn_std = std(sdnn_vals);
sdnn_ci = t_score * (sdnn_std / sqrt(n));
fprintf("\n--- SDNN ---\nMean = %.2f ms | 95%% CI = [%.2f, %.2f] ms\n", ...
        sdnn_mean, sdnn_mean - sdnn_ci, sdnn_mean + sdnn_ci);

% RMSSD CI
rmssd_mean = mean(rmssd_vals);
rmssd_std = std(rmssd_vals);
rmssd_ci = t_score * (rmssd_std / sqrt(n));
fprintf("--- RMSSD ---\nMean = %.2f ms | 95%% CI = [%.2f, %.2f] ms\n", ...
        rmssd_mean, rmssd_mean - rmssd_ci, rmssd_mean + rmssd_ci);
