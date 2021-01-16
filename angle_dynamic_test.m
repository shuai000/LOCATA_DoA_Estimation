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
% index_all = [1 2; 5 6; 5 7; 6 7; 4 5; 11 12];
index_all = [1 2; 5 6; 5 7; 6 7; 13 15; 14 15];

[true_delay, timestamps, audio_array, source, mic_positions, target, h, R, azimuth_truth, elevation_truth] = main2(array_dir, this_array, index_all);

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

%% Location estimation
Ne = length(index_array);
azimuth_estimation = zeros(1, Ne);
elevation_estimation = zeros(1, Ne);
position_estimation = zeros(3, Ne);
r = zeros(1, Ne);
for i=1:Ne
    [azimuth_estimation(i), elevation_estimation(i), r(i)] = find_angle(h, R, target.talker5.position(:, i));
end
figure;
subplot 211
plot(rad2deg(azimuth_estimation)); hold on;
plot(rad2deg(azimuth_truth));
legend('estimation', 'truth');

subplot 212
plot(rad2deg(elevation_estimation)); hold on;
plot(rad2deg(elevation_truth));
legend('estimation', 'truth');
% h_est_az = plot(results.source(src_idx).timestamps, rad2deg(wrapToPi(results.source(src_idx).azimuth)), 'bo');
% % Plot the ground truth
% h_true = plot(in.timestamps, rad2deg(truth.source.(src).azimuth), 'rx');
