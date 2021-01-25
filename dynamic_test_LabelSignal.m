%%% Localization
%%% 06/01/2021
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

% % addpath C:\Users\61450\OneDrive - RMIT University\Documents\Parametric_LAP -- MicrophoneArray';
% audio_data = importdata('D:\LOCATA\dev\task3\recording1\dicit\audio_array_dicit.wav');
% recording2, talker2
% recording1, talker5
% recording3, talker1

% Select the array sensor index
m1 = 5;
m2 = 6;
m3 = 7;
m4 = 9;

% index_all = [m1 m2; m1 m3; m1 m4; m2 m3; m2 m4; m3 m4];
% index_all = [1 2; 5 6; 5 7; 6 7; 4 5; 11 12];
index_all = [1 2; 5 6; 5 7; 6 7; 13 15; 14 15];

[true_delay, timestamps, audio_array, source, mic_positions, target, h, R, azimuth_truth, elevation_truth] = main2(array_dir, this_array, index_all, field);

delay_max = max(max(true_delay));

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

%% This is to contatecate the segment data
% Na = sum(label);
% label_data = zeros(1, Na);
% label_index = zeros(1, Na);
% counter = 0;
% for i=1:length(source_audio)
%     if label(i) == 1
%         counter = counter + 1;
%         label_data(counter) = source_audio(i);
%         label_index(counter) = i;
%     end
% end
%
% t_label = t(label_index);  % the corresponding time for the labelled activity signal
% data = data(label_index, :); % extract the labelled data
%
% figure;
% plot(t_label, label_data, 'r.'); hold on;
% plot(t, source_audio, 'b');
% legend('labelled signal', 'original signal');
% xlabel('time, t (seconds)');

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

%% Parameters for delay estimation:
scale = 6;  % test scales,
window = 1000;  % large window size? during window constant delay, 4000
orders = 2:10;  % test order, which one it picks, remove 1

N_segment = length(data_segment);
N_pair = size(index_all, 1);
delay_estimate = cell(1, N_segment);
order = cell(1, N_segment);

for ks=1:N_segment
    
    delay_estimate{ks} = zeros(size(data_segment(ks).data, 1), N_pair);
    order{ks} = zeros(1, N_pair);
    
    for i=1:N_pair
        index_tx = index_all(i, 1);
        index_rx = index_all(i, 2);
        signal_tx = data_segment(ks).data(:, index_tx)';
        signal_rx = data_segment(ks).data(:, index_rx)';
        
        %% Load transmitted and received signal:
        [delayEst,Order,~] = MultiScale_LAP_Param(signal_rx, signal_tx, scale,window,orders);
        order{ks}(i) = Order;
        signalEst = imshift(signal_tx,1i.*delayEst);
        MSE = mean(abs((signal_rx - signalEst)).^2)/mean(signal_rx.^2);
        
        fprintf('\t\t\t\tNo Iterations:\t MSE = %1.5f, Delay Difference = %1.4f\n', MSE, mean(abs(diff(delayEst)))/4);
        
        delay_estimate{ks}(:, i) =  delayEst';
        
    end
end

% error_o = abs(delay_estimate(index_array, :)' - true_delay(:, true_array));
% error = reshape(error_o, 1, size(error_o, 1) * size(error_o, 2));
% fprintf('\t MaxError = %1.5f, Mean Error = %1.4f\n', max(error), mean(mean(error)));

% figure;
% % plot error
% for i=1:size(index_all, 1)
%     plot(time_array_estimation, abs(delay_estimate(index_array,i)' - true_delay(i, true_array)));
%     hold on;
% end
% legend('1', '2', '3', '4', '5', '6');
% xlabel('time (seconds)');
% ylabel('absolute estimation error (sample)');

figure;
subplot 611
plot(timestamps, true_delay(1, :), 'k'); hold on;
for k=1:N_segment
    plot(data_segment(k).time, delay_estimate{k}(:, 1), 'r'); hold on;
end
ylim([-50, 50]);
legend('truth', 'estimation');

subplot 612
plot(timestamps, true_delay(2, :), 'k'); hold on;
for k=1:N_segment
    plot(data_segment(k).time, delay_estimate{k}(:, 2), 'r'); hold on;
end
ylim([-50, 50]);

subplot 613
plot(timestamps, true_delay(3, :), 'k'); hold on;
for k=1:N_segment
    plot(data_segment(k).time, delay_estimate{k}(:, 3), 'r'); hold on;
end
ylim([-50, 50]);
ylabel('delay (sample)');

subplot 614
plot(timestamps, true_delay(4, :), 'k'); hold on;
for k=1:N_segment
    plot(data_segment(k).time, delay_estimate{k}(:, 4), 'r'); hold on;
end
ylim([-50, 50]);
subplot 615
plot(timestamps, true_delay(5, :), 'k'); hold on;
for k=1:N_segment
    plot(data_segment(k).time, delay_estimate{k}(:, 5), 'r'); hold on;
end
ylim([-50, 50]);
subplot 616
plot(timestamps, true_delay(6, :), 'k'); hold on;
for k=1:N_segment
    plot(data_segment(k).time, delay_estimate{k}(:, 6), 'r'); hold on;
end
ylim([-50, 50]);
xlabel('Time, $t$, [s]', 'interpreter', 'latex');

% %% Location estimation
azimuth_estimation = cell(1, N_segment);
elevation_estimation = cell(1, N_segment);
position_estimation = cell(1, N_segment);
r = cell(1, N_segment);

for k=1:N_segment
    k
    Ne = length(meta_segment(k).local_index);
    azimuth_estimation{k} = zeros(1, Ne);
    elevation_estimation{k} = zeros(1, Ne);
    position_estimation{k} = zeros(3, Ne);
    r{k} = zeros(1, Ne);
    
    current_azimuth = zeros(1, Ne);
    current_elevation = zeros(1, Ne);
    current_r = zeros(1, Ne);
    
    for i=1:Ne
        estimated_delay = delay_estimate{k}(meta_segment(k).local_index(i), :);
        % estimated_delay = true_delay(:, i)';
        [x, y, z, error_min] = location_search2(estimated_delay, mic_positions, index_all);
        position_estimation{k}(:, i) = [x, y, z]';
        [current_azimuth(i), current_elevation(i), current_r(i)] = find_angle(h, R, position_estimation{k}(:, i));
    end
    azimuth_estimation{k} = current_azimuth;
    elevation_estimation{k} = current_elevation;
    r{k} = current_r;
end

mean_error = zeros(1, N_segment);
for k=1:N_segment
    mean_error(k) = mean(abs(rad2deg(azimuth_estimation{k}) - rad2deg(azimuth_truth(meta_segment(k).true_index))));
end

fprintf('\t\t Mean Error = %1.5f \n', mean(mean_error));

figure;
subplot 211
plot(time_array_estimation, rad2deg(azimuth_truth(true_array)), 'rx'); hold on;
for k=1:N_segment
    plot(data_segment(k).time(meta_segment(k).local_index), rad2deg(azimuth_estimation{k}), 'bo'); hold on;
end
legend('ground truth', 'estimate', 'interpreter', 'latex');
ylabel('Azimuth [Deg]', 'interpreter', 'latex');
ylim([-180, 180]);
grid on;

subplot 212
plot(time_array_estimation, rad2deg(elevation_truth(true_array)), 'rx'); hold on;
for k=1:N_segment
    plot(data_segment(k).time(meta_segment(k).local_index), rad2deg(elevation_estimation{k}), 'bo'); hold on;
end

legend('estimate', 'ground truth', 'interpreter', 'latex');
ylabel('Elevation [Deg]', 'interpreter', 'latex');
ylim([0, 180]);
xlabel('Time, t (seconds)', 'interpreter', 'latex');
grid on;

figure;
subplot 311
plot(time_array_estimation,target.(field).position(1, true_array)); hold on;
for k=1:N_segment
    plot(data_segment(k).time(meta_segment(k).local_index), position_estimation{k}(1, :)); hold on;
end
legend('truth', 'estimation', 'interpreter', 'latex');
ylabel('x (m)', 'interpreter', 'latex');
ylim([-1.5, 1.5]);

subplot 312
plot(time_array_estimation,target.(field).position(2, true_array)); hold on;
for k=1:N_segment
    plot(data_segment(k).time(meta_segment(k).local_index), position_estimation{k}(2, :)); hold on;
end
legend('truth', 'estimation',  'interpreter', 'latex');
ylabel('y (m)', 'interpreter', 'latex');
ylim([-2, 2]);

subplot 313
plot(time_array_estimation,target.(field).position(3, true_array)); hold on;
for k=1:N_segment
    plot(data_segment(k).time(meta_segment(k).local_index), position_estimation{k}(3, :)); hold on;
end
legend('truth', 'estimation', 'interpreter', 'latex');
ylabel('z (m)', 'interpreter', 'latex');
ylim([0, 1.6]);
