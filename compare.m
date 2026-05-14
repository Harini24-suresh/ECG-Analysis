% Benchmark Script: Compare Pan-Tompkins with Ground Truth
records = {'a01', 'a02', 'a04','a19','b04','c04'};
fs = 100;

for i = 1:length(records)
    fprintf("\n=== Benchmarking: %s ===\n", records{i});
    benchmark_rpeak_detection(records{i}, fs);
end
