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

load delay_estimate.mat

% Select the array sensor index
m1 = 5;
m2 = 6;
m3 = 7;
m4 = 9;

% index_all = [m1 m2; m1 m3; m1 m4; m2 m3; m2 m4; m3 m4];
index_all = [1 2; 5 6; 5 7; 6 7; 4 5; 11 12];

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

k=1;
delay_test = delay_estimate(1, :);
[x, y, z, min_error] = location_search2(delay_test, mic_positions, index_all);

target.talker5.position(:, k)
[x, y, z]
min_error

% figure;
% plot(index_array, 'r-o');

% index = [];
% for i = 1:length(timestamps)
%     current_index = find(t == timestamps(i));
%     index = [index, current_index];
% end

% sound(source_audio, fs);



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
