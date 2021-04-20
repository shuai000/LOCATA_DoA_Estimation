%%% Localization
%%% 06/01/2021
%%% Shuai SUN
clear all;
close all;
clc;

array_dir = 'D:\LOCATA\dev\task3\recording3\dicit';
this_array = 'dicit';
field = 'talker1';
% % addpath C:\Users\61450\OneDrive - RMIT University\Documents\Parametric_LAP -- MicrophoneArray';
% audio_data = importdata('D:\LOCATA\dev\task3\recording1\dicit\audio_array_dicit.wav');
% recording2, talker2
% recording1, talker5
% recording3, talker1

index_all = [1 2; 5 6; 5 7; 6 7; 13 15; 14 15];

[true_delay, timestamps, audio_array, source, mic_positions, target, h, R, azimuth_truth, elevation_truth] = main2(array_dir, this_array, index_all, field);

delay_max = max(max(true_delay));

source_audio = source.data.(field);

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

sound(source_audio, fs);

%% Parameters for delay estimation:
scale = 6;  % test scales,
window = 10000;  % large window size? during window constant delay 
orders = 4;
order_vector = [3, 4, 4, 4, 6, 5];

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

error_o = abs(delay_estimate(index_array, :)' - true_delay(:, true_array));
error = reshape(error_o, 1, size(error_o, 1) * size(error_o, 2));
fprintf('\t MaxError = %1.5f, Mean Error = %1.4f\n', max(error), mean(mean(error)));

figure; 
% plot error
for i=1:size(index_all, 1)
    plot(time_array_estimation, abs(delay_estimate(index_array,i)' - true_delay(i, true_array)));
    hold on;
end
legend('1', '2', '3', '4', '5', '6');
xlabel('time (seconds)');
ylabel('absolute estimation error (sample)');

figure;
subplot 611
plot(t(valid_idx), delay_estimate(valid_idx, 1), 'r'); hold on;
plot(timestamps, true_delay(1, :), 'k');
legend('estimation', 'truth');

subplot 612
plot(t(valid_idx), delay_estimate(valid_idx, 2), 'r'); hold on;
plot(timestamps, true_delay(2, :), 'k');

subplot 613
plot(t(valid_idx), delay_estimate(valid_idx, 3), 'r'); hold on;
plot(timestamps, true_delay(3, :), 'k');
ylabel('delay (sample)');

subplot 614
plot(t(valid_idx), delay_estimate(valid_idx, 4), 'r'); hold on;
plot(timestamps, true_delay(4, :), 'k');

subplot 615
plot(t(valid_idx), delay_estimate(valid_idx, 5), 'r'); hold on;
plot(timestamps, true_delay(5, :), 'k');

subplot 616
plot(t(valid_idx), delay_estimate(valid_idx, 6), 'r'); hold on;
plot(timestamps, true_delay(6, :), 'k');
xlabel('Time, $t$, [s]', 'interpreter', 'latex');

% %% Location estimation
Ne = length(index_array);
azimuth_estimation = zeros(1, Ne);
elevation_estimation = zeros(1, Ne);
position_estimation = zeros(3, Ne);
r = zeros(1, Ne);
parfor i=1:Ne
    i
    estimated_delay = delay_estimate(index_array(i), :);
%     estimated_delay = true_delay(:, i)';
    [x, y, z, error_min] = location_search2(estimated_delay, mic_positions, index_all);
    position_estimation(:, i) = [x, y, z]';
    [azimuth_estimation(i), elevation_estimation(i), r(i)] = find_angle(h, R, position_estimation(:, i));
end
figure;
subplot 211
plot(time_array_estimation, rad2deg(azimuth_estimation), 'bo'); hold on;
plot(time_array_estimation, rad2deg(azimuth_truth(true_array)), 'rx');
legend('estimate', 'ground truth', 'interpreter', 'latex');
ylabel('Azimuth [Deg]', 'interpreter', 'latex');
ylim([-180, 180]);
grid on;

subplot 212
plot(time_array_estimation, rad2deg(elevation_estimation), 'bo'); hold on;
plot(time_array_estimation, rad2deg(elevation_truth(true_array)), 'rx');
legend('estimate', 'ground truth', 'interpreter', 'latex');
ylabel('Elevation [Deg]', 'interpreter', 'latex');
ylim([0, 180]);
xlabel('Time, t (seconds)', 'interpreter', 'latex');
grid on;

figure;
subplot 311
plot(time_array_estimation, position_estimation(1, :)); hold on;
plot(time_array_estimation,target.(field).position(1, true_array));
legend('estimation', 'truth', 'interpreter', 'latex');
ylabel('x (m)', 'interpreter', 'latex');
ylim([-1.5, 1.5]);

subplot 312
plot(time_array_estimation, position_estimation(2, :)); hold on;
plot(time_array_estimation,target.(field).position(2, true_array));
legend('estimation', 'truth', 'interpreter', 'latex');
ylabel('y (m)', 'interpreter', 'latex');
ylim([-2, 2]);

subplot 313
plot(time_array_estimation, position_estimation(3, :)); hold on;
plot(time_array_estimation,target.(field).position(3, true_array));
legend('estimation', 'truth', 'interpreter', 'latex');
ylabel('z (m)', 'interpreter', 'latex');
ylim([0, 1.6]);
