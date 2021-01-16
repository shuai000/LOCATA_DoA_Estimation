%%% Localization
%%% 06/01/2021
%%% Shuai SUN
clear all;
close all;
clc;

array_dir = 'D:\LOCATA\dev\task3\recording1\dicit';
this_array = 'dicit';
% % addpath C:\Users\61450\OneDrive - RMIT University\Documents\Parametric_LAP -- MicrophoneArray';
% audio_data = importdata('D:\LOCATA\dev\task3\recording1\dicit\audio_array_dicit.wav');

% Select the array sensor index
m1 = 5;
m2 = 6;
m3 = 7;
m4 = 9;

% index_all = [m1 m2; m1 m3; m1 m4; m2 m3; m2 m4; m3 m4];
index_all = [12 15; 14 15; 13 15; 11 15; 1 8; 2 8];

[true_delay, timestamps, audio_array, source, mic_positions, target] = main2(array_dir, this_array, index_all);

source_audio = source.data.talker5;

data = audio_array.data.dicit;
fs = audio_array.fs;

figure;
t = linspace(0,(size(data,1)-1)/fs, size(data,1));
valid_idx = find(t <= timestamps(end));  % data that is longer than the opitracker is dropped
plot(t(valid_idx), source_audio(valid_idx));
xlabel('Time, $t$, [s]', 'interpreter', 'latex');
ylabel('Amplitude', 'interpreter', 'latex');
grid on;

Ns = length(timestamps);
index_array = [];
true_array = [];
for i = 1:Ns
    % find index to extract delay estimation, not perfertly synchronized
    current_index = find_index(timestamps(i), t);
    if current_index > 0
        index_array = [index_array, current_index];
        true_array = [true_array, i];
    end
end
time_array_estimation = t(index_array);

% k=1;
% [x, y, z] = location_search(true_delay(:, k), mic_positions, index_all);
% target.talker5.position(:, k)
% figure;
% plot(index_array, 'r-o');

% index = [];
% for i = 1:length(timestamps)
%     current_index = find(t == timestamps(i));
%     index = [index, current_index];
% end

% sound(source_audio, fs);

%% Parameters for delay estimation:
scale = 8;
window = 1000;
orders = 4;
order_vector = [6, 4, 4, 4, 4, 6];

delay_estimate = zeros(size(data,1), size(index_array, 1));
order = [];

for i=1:size(index_all, 1)
    index_tx = index_all(i, 1);
    index_rx = index_all(i, 2);
    signal_tx = data(:, index_tx)';
    signal_rx = data(:, index_rx)';
    
    %% Load transmitted and received signal:
    [delayEst,Order,~] = MultiScale_LAP_Param(signal_rx, signal_tx, scale,window,order_vector(i));
    order = [order, Order];
    signalEst = imshift(signal_tx,1i.*delayEst);
    MSE = mean(abs((signal_rx - signalEst)).^2)/mean(signal_rx.^2);
    
    fprintf('\t\t\t\tNo Iterations:\t MSE = %1.5f, Delay Difference = %1.4f\n', MSE, mean(abs(diff(delayEst)))/4);
    
    delay_estimate(:, i) =  delayEst';
    
    %     figure;
    %     plot(t(valid_idx), delay_estimate(valid_idx, i), 'r'); hold on;
    %     plot(timestamps, true_delay(i, :), 'k.');
    
    %     figure;
    %     plot(signal_tx, 'r'); hold on;
    %     plot(signal_rx, 'k'); legend("s1", "s2");
end

figure; 
% plot error
for i=1:size(index_all, 1)
    plot(time_array_estimation, abs(delay_estimate(index_array,i)' - true_delay(i, true_array)));
    hold on;
end

legend('1', '14-15', '13-15', '11-15', '1-8', '2-8');
xlabel('time (seconds)');
ylabel('absolute estimation error (sample)');

figure;
subplot 611
plot(t(valid_idx), delay_estimate(valid_idx, 1), 'r'); hold on;
plot(timestamps, true_delay(1, :), 'k');
legend('estimation', 'truth');
title('Array: 12 - 15');

subplot 612
plot(t(valid_idx), delay_estimate(valid_idx, 2), 'r'); hold on;
plot(timestamps, true_delay(2, :), 'k');
title('Array: 14 -15');

subplot 613
plot(t(valid_idx), delay_estimate(valid_idx, 3), 'r'); hold on;
plot(timestamps, true_delay(3, :), 'k');
ylabel('delay (sample)');
title('Array: 13-15');

subplot 614
plot(t(valid_idx), delay_estimate(valid_idx, 4), 'r'); hold on;
plot(timestamps, true_delay(4, :), 'k');
title('Array: 11 - 15');

subplot 615
plot(t(valid_idx), delay_estimate(valid_idx, 5), 'r'); hold on;
plot(timestamps, true_delay(5, :), 'k');
title('Array: 1 - 8');

subplot 616
plot(t(valid_idx), delay_estimate(valid_idx, 6), 'r'); hold on;
plot(timestamps, true_delay(6, :), 'k');
xlabel('Time, $t$, [s]', 'interpreter', 'latex');
title('Array: 2 - 8');

%% Location estimation
azimuth_estimation = zeros(1, Ns);
elevation_estimation = zeros(1, Ns);
position_estimation = zeros(3, Ns);

% for i=1
%     estimated_delay = delay_estimate(index_array(i), :);
%     [x, y, z] = location_search(estimated_delay, mic_positions, index_all);
% end

% h_est_az = plot(results.source(src_idx).timestamps, rad2deg(wrapToPi(results.source(src_idx).azimuth)), 'bo');
% % Plot the ground truth
% h_true = plot(in.timestamps, rad2deg(truth.source.(src).azimuth), 'rx');
