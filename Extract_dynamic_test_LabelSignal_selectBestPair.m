%%% Script for audio soucrce Localization
%%% TDoA approach (time delay is estimated by Parametric LAP algorithm)
%%% segment-based apporach, modified on 27/01/2021
%%% Shuai SUN
clear all;
close all;
clc;

index_all = [7 12];


this_array = 'dicit';

array_dirc = cell(1, 3);
array_dirc{1} = 'D:\LOCATA\dev\task3\recording1\dicit';
array_dirc{2} = 'D:\LOCATA\dev\task3\recording2\dicit';
array_dirc{3} = 'D:\LOCATA\dev\task3\recording3\dicit';

fieldc = cell(1, 3);
fieldc{1} = 'talker5';
fieldc{2} = 'talker2';
fieldc{3} = 'talker1';

label_dirc = cell(1, 3);
label_dirc{1} = 'D:\LOCATA\dev\task3\recording1\dicit\VAD_dicit_talker5.txt';
label_dirc{2} = 'D:\LOCATA\dev\task3\recording2\dicit\VAD_dicit_talker2.txt';
label_dirc{3} = 'D:\LOCATA\dev\task3\recording3\dicit\VAD_dicit_talker1.txt';

% recording2, talker2
% recording1, talker5
% recording3, talker1

for mi=3
    field = fieldc{mi};
    label_dir = label_dirc{mi};
    array_dir = array_dirc{mi};
    
    file_t = importdata(label_dir);
    label = file_t.data;
    
    [true_delay, timestamps, audio_array, source, mic_positions, target, h, R, azimuth_truth, elevation_truth] = main2(array_dir, this_array, index_all, field);
    
    delay_max = max(max(true_delay)); % compute the possible true maximum delay
    
    source_audio = source.data.(field);
    
    data = audio_array.data.dicit;
    fs = audio_array.fs;
    
    t = linspace(0,(size(data,1)-1)/fs, size(data,1));
    
    t_step = (size(data,1)-1)/fs;
    
    [data_segment, meta_segment] = find_segment(data, label, t, timestamps);
    
    figure;
    plot(data_segment(1).time, data_segment(1).data(:, 7)); hold on;
    plot(data_segment(1).time, data_segment(1).data(:, 12));
    
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
    
    %% Parameters for delay estimation:
    scale = 6;  % test scales,
    window = 4000;  % large window size? during window constant delay, 4000
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
            
            %             fprintf('\t\t\t\tNo Iterations:\t MSE = %1.5f, Delay Difference = %1.4f\n', MSE, mean(abs(diff(delayEst)))/4);
            
            delay_estimate{ks}(:, i) =  delayEst';
            
        end
    end
    
    figure;
    yyaxis left;
    plot(timestamps, true_delay(1, :), 'k'); hold on;
    for k=1:N_segment
        plot(data_segment(k).time, delay_estimate{k}(:, 1), 'r'); hold on;
    end
    ylim([-50, 50]);
    yyaxis right;
    plot(data_segment(1).time, data_segment(1).data(:, 7)); hold on;
    plot(data_segment(1).time, data_segment(1).data(:, 12));
    
    figure;
    t1 = 6900;
    t2 = t1 + 600;
    plot(data_segment(1).time(t1:t2), data_segment(1).data(t1:t2, 7), 'r', 'Markersize', 3); hold on;
    plot(data_segment(1).time(t1:t2), data_segment(1).data(t1:t2, 12), 'k', 'Markersize', 3);
  
    
    figure;
    t1 = 6900;
    t2 = t1 + 600;
    dt = 8;
    plot(data_segment(1).time(t1:t2), data_segment(1).data(t1+dt:t2+dt, 7), 'r', 'Markersize', 3); hold on;
    plot(data_segment(1).time(t1:t2), data_segment(1).data(t1:t2, 12), 'k', 'Markersize', 3);
  
  
    figure;
    t1 = 6000;
    t2 = t1 + 6000;
    plot(data_segment(1).time(t1:t2), data_segment(1).data(t1:t2, 7), 'r', 'Markersize', 3); hold on;
    plot(data_segment(1).time(t1:t2), data_segment(1).data(t1:t2, 12), 'k', 'Markersize', 3);
    xlabel('Time (seconds)', 'interpreter', 'latex');
    xlabel('Amplitude', 'interpreter', 'latex');
    legend('signal1', 'signal2', 'interpreter', 'latex');
    ylim([-.05, 0.08]);
    
    figure;
    t1 = 6000;
    t2 = t1 + 6000;
    dt = 8;
    plot(data_segment(1).time(t1:t2), data_segment(1).data(t1+dt:t2+dt, 7), 'r', 'Markersize', 3); hold on;
    plot(data_segment(1).time(t1:t2), data_segment(1).data(t1:t2, 12), 'k', 'Markersize', 3);
   xlabel('Time (seconds)', 'interpreter', 'latex');
    xlabel('Amplitude', 'interpreter', 'latex');
    legend('signal1', 'signal2', 'interpreter', 'latex');
      ylim([-.05, 0.08]);
    
    % %% Location estimation
    azimuth_estimation = cell(1, N_segment);
    elevation_estimation = cell(1, N_segment);
    position_estimation = cell(1, N_segment);
    r = cell(1, N_segment);
    
    for k=1:N_segment
        Ne = length(meta_segment(k).local_index);
        azimuth_estimation{k} = zeros(1, Ne);
        elevation_estimation{k} = zeros(1, Ne);
        position_estimation{k} = zeros(3, Ne);
        r{k} = zeros(1, Ne);
        
        current_azimuth = zeros(1, Ne);
        current_elevation = zeros(1, Ne);
        current_r = zeros(1, Ne);
        
%         for i=1:Ne
%             estimated_delay = delay_estimate{k}(meta_segment(k).local_index(i), :);
%             %             estimated_delay = true_delay(:, meta_segment(k).true_index(i))';
%             [x, y, z, error_min] = location_search2(estimated_delay, mic_positions, index_all);
%             position_estimation{k}(:, i) = [x, y, z]';
%             [current_azimuth(i), current_elevation(i), current_r(i)] = find_angle(h, R, position_estimation{k}(:, i));
%         end
        azimuth_estimation{k} = current_azimuth;
        elevation_estimation{k} = current_elevation;
        r{k} = current_r;
    end
    
    mean_error = zeros(1, N_segment);
    for k=1:N_segment
        mean_error(k) = mean(abs(rad2deg(azimuth_estimation{k}) - rad2deg(azimuth_truth(meta_segment(k).true_index))));
    end
    
    fprintf('\t\t Mean Error for this recording = %1.5f \n', mean(mean_error));
    
    %% Extract data 21/02/2021
    ke = 1;
    Ne = length(meta_segment(ke).local_index);
    
    truth_delay = zeros(N_pair, Ne);
    truth_angle = zeros(2, Ne);
    estimate_delay = zeros(N_pair, Ne);
    estimate_angle = zeros(2, Ne);
    
    truth_angle(1, :) = rad2deg(azimuth_truth(meta_segment(ke).true_index));
    truth_angle(2, :) = rad2deg(elevation_truth(meta_segment(ke).true_index));
    
    for i=1:N_pair
        truth_delay(i, :) = true_delay(i, meta_segment(ke).true_index)';
        estimate_delay(i, :) = delay_estimate{ke}(meta_segment(ke).local_index, i)';
    end
    
    estimate_angle(1, :) =  rad2deg(azimuth_estimation{ke});
    estimate_angle(2, :) = rad2deg(elevation_estimation{ke});
    
    save('D:\LOCATA\sap_locata_io\space32cm.mat', 'truth_delay', 'truth_angle', 'estimate_delay', 'estimate_angle');
    
end


%%
% load('center7.mat');
%
% % sampling frequency:
% fs = 48000;
%
% % distance between mics:
% d = [0.64,0.32,0.16,0.08,0.04, 0.04,0.08, 0.16, 0.32,0.64].';
%
% % speed of sound:
% c = 343;
%
% % estimation of azimuth:
% azimuth_est = zeros(1,205);
% for index = 1:205
%     tau = estimate_delay(:,index)./fs;
%     azimuth_est(index) = asin((d\tau)*c)*180/pi;
% end