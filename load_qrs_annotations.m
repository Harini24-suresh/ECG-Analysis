% benchmark_runner.m
fs = 100;  % sampling frequency
records = {'a02', 'a04', 'a19'};  % records for which you have .qrs files

for i = 1:length(records)
    benchmark_rpeak_detection(records{i}, fs);
end
