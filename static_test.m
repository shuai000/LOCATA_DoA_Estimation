%%% Localization
%%% 06/01/2021
%%% Shuai SUN
clear all;
close all;
clc;

array_dir = 'D:\LOCATA\dev\task1\recording3\dicit';
this_array = 'dicit';
% addpath C:\Users\61450\OneDrive - RMIT University\Documents\Parametric_LAP -- MicrophoneArray';
audio_data = importdata('D:\LOCATA\dev\task1\recording3\dicit\audio_array_dicit.wav');

data = audio_data.data;
fs = audio_data.fs;

% Select the array sensor index
m1 = 5;
m2 = 6;
m3 = 7;
m4 = 9;

% index_all = [m1 m2; m1 m3; m1 m4; m2 m3; m2 m4; m3 m4];
index_all = [1 2; 5 6; 5 7; 6 7; 4 5; 11 12];

true_delay = main1(array_dir, this_array, index_all);


%% Parameters for delay estimation:
scale = 8;
window = 1000;
orders = 1:10;

delay_estimate = zeros(1, size(index_all, 1));
order = [];

for i=1:size(index_all, 1)
    index_tx = index_all(i, 1);
    index_rx = index_all(i, 2);
    signal_tx = data(:, index_tx)';
    signal_rx = data(:, index_rx)';
    
    %% Load transmitted and received signal:
    [delayEst,Order,~] = MultiScale_LAP_Param(signal_rx, signal_tx, scale,window,orders);
    order = [order, Order];
    signalEst = imshift(signal_tx,1i.*delayEst);
    MSE = mean(abs((signal_rx - signalEst)).^2)/mean(signal_rx.^2);
    
    fprintf('\t\t\t\tNo Iterations:\t MSE = %1.5f, Delay Difference = %1.4f\n', MSE, mean(abs(diff(delayEst)))/4);
    
    delay_estimate(i) =  delayEst(end);
    
    % figure;
    % plot(signal_tx, 'r'); hold on;
    % plot(signal_rx, 'k'); legend("s1", "s2");
end

delay_estimate - true_delay
