%%% Script for extracting required angles for Chris
%%% TDoA approach (time delay is estimated by Parametric LAP algorithm)
%%% segment-based apporach, modified on 18/02/2021
%%% *** TEST for delay estimation (which pairs are best)
%%% Shuai SUN
clear all;
close all;
clc;

array_dir = 'D:\LOCATA\dev\task3\recording3\dicit';
this_array = 'dicit';
field = 'talker1';

label_dir = 'D:\LOCATA\dev\task3\recording3\dicit\VAD_dicit_talker1.txt';
file_t = importdata(label_dir);
label = file_t.data;

% recording2, talker2
% recording1, talker5
% recording3, talker1

% index_all = [1 2; 5 6; 5 7; 6 7; 13 15; 14 15];  % 7 9 error 15
% 32cm
% index_all = [1 8; 1 2; 2 3; 3 7; 7 12; 12 13; 13 14; 14 15];
% 16cm
% index_all = [3 4; 4 7; 5 9; 6 10; 7 11; 11 12];
% 8cm
% index_all = [4 5; 5 7; 6 9; 7 10; 10 11];
% 4cm
% index_all = [5 6; 6 7; 7 9; 9 10];
% 64cm
% index_all = [1 3; 2 7; 3 12; 7 13; 12 14; 13 15; 2 8];
% centered on 7:
% index_all = [4 7; 5 7; 6 7; 7 9; 7 10; 7 11];

index_all = [2 7; 3 7; 4 7; 5 7; 6 7; 7 9; 7 10; 7 11; 7 12; 7 13];

% index_all = [1 2; 2 3; 3 7; 7 12; 12 13; 13 14];


[true_delay, timestamps, audio_array, source, mic_positions, target, h, R, azimuth_truth, elevation_truth] = main2(array_dir, this_array, index_all, field);

delay_max = max(max(true_delay)); % compute the possible true maximum delay

source_audio = source.data.(field);

data = audio_array.data.dicit;
fs = audio_array.fs;

t = linspace(0,(size(data,1)-1)/fs, size(data,1));

[data_segment, meta_segment] = find_segment(data, label, t, timestamps);

figure;
valid_idx = find(t <= timestamps(end));  % data that is longer than the opitracker is dropped
plot(t(valid_idx), source_audio(valid_idx)); hold on;
for k=1:length(data_segment)
    k1 = data_segment(k).index_interval(1);
    k2 = data_segment(k).index_interval(2);
    plot(data_segment(k).time, source_audio(k1:k2), 'r.'); hold on;
end
xlabel('Time, $t$, [s]', 'interpreter', 'latex');
ylabel('Amplitude', 'interpreter', 'latex');
grid on;

figure;
subplot 211
valid_idx = find(t <= timestamps(end));  % data that is longer than the opitracker is dropped
plot(t(valid_idx), source_audio(valid_idx)); hold on;
for k=1:length(data_segment)
    k1 = data_segment(k).index_interval(1);
    k2 = data_segment(k).index_interval(2);
    plot(data_segment(k).time, source_audio(k1:k2), 'r.'); hold on;
end
xlabel('Time, $t$, [s]', 'interpreter', 'latex');
ylabel('Amplitude', 'interpreter', 'latex');
grid on;

subplot 212
for k=1:length(meta_segment)
    current_time = timestamps(meta_segment(k).true_index);
    data_toplot = data_segment(k).data(meta_segment(k).local_index, 1);
    plot(current_time, data_toplot, 'm'); hold on;
end
xlabel('Time, $t$, [s]', 'interpreter', 'latex');
ylabel('Ground Truth time', 'interpreter', 'latex');

Ns = length(timestamps);
index_array = [];
true_array = [];
for i = 1:Ns
    % find index to extract delay estimation, not perfertly synchronized
    current_index = find_index(timestamps(i), t);
    if current_index > 0
        index_array = [index_array, current_index]; % w.r.t. global data index, (long)
        true_array = [true_array, i]; % w.r.t. timestamp (short)
    end
end
time_array_estimation = t(index_array);

% sound(source_audio, fs);

N_segment = length(data_segment);
N_pair = size(index_all, 1);
delay_estimate = cell(1, N_segment);
order = cell(1, N_segment);

k = 1; % the first segment
Ne = length(meta_segment(k).local_index);
true_delay_pair = zeros(N_pair, Ne);
for i=1:N_pair
    true_delay_pair(i, :) = true_delay(i, meta_segment(k).true_index)';
end

true_angle_array = zeros(2, Ne);
true_angle_array(1, :) = rad2deg(azimuth_truth(meta_segment(k).true_index));
true_angle_array(2, :) = rad2deg(elevation_truth(meta_segment(k).true_index));

%% estimation
scale = 6;  % test scales,
window = 4000;  % large window size? during window constant delay, 4000
orders = 2:10;  % test order, which one it picks, remove 1

estimation_delay = zeros(N_pair, Ne);

for i=1:N_pair
    index_tx = index_all(i, 1);
    index_rx = index_all(i, 2);
    signal_tx = data_segment(k).data(:, index_tx)';
    signal_rx = data_segment(k).data(:, index_rx)';
    
    %% Load transmitted and received signal:
    [delayEst,Order,~] = MultiScale_LAP_Param(signal_rx, signal_tx, scale,window,orders);
    signalEst = imshift(signal_tx,1i.*delayEst);
    MSE = mean(abs((signal_rx - signalEst)).^2)/mean(signal_rx.^2);
    
    estimation_delay(i, :) =  delayEst(meta_segment(k).local_index);
    
end

error_tolook = estimation_delay - true_delay_pair;

%% Azimuth estimation
estimation_angle = zeros(2, Ne);
for i=1:Ne
     estimated_delay = estimation_delay(:, i)';
    [x, y, z, error_min] = location_search2(estimated_delay, mic_positions, index_all);
    [estimation_angle(1, i), estimation_angle(2, i), ~] = find_angle(h, R, [x, y, z]');
end

save('D:\LOCATA\sap_locata_io\center7.mat', 'true_angle_array', 'true_delay_pair', 'estimation_delay', 'estimation_angle');

figure;

for i=1:N_pair
    plot(estimation_delay(i, :) - true_delay_pair(i, :)); hold on;
end

figure;
subplot 211
plot(time_array_estimation, rad2deg(azimuth_truth(true_array)), 'bx'); hold on;
plot(data_segment(k).time(meta_segment(k).local_index), rad2deg(azimuth_truth(meta_segment(k).true_index)), 'rx');
ylabel('Azimuth');

subplot 212
plot(time_array_estimation, rad2deg(elevation_truth(true_array)), 'bx'); hold on;
plot(data_segment(k).time(meta_segment(k).local_index), rad2deg(elevation_truth(meta_segment(k).true_index)), 'rx');
ylabel('Elevation');
xlabel('Time (seconds)');

figure;
subplot 311
plot(time_array_estimation,target.(field).position(1, true_array)); hold on;

subplot 312
plot(time_array_estimation,target.(field).position(2, true_array)); hold on;

subplot 313
plot(time_array_estimation,target.(field).position(3, true_array));

