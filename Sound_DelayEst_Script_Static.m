% Script for estimating the delay between microphone sound
% The delay is estimated assuming a parametric structure.
% 20/11/2020 Shuai SUN

clc; clear all; close all;

load signal_tx.mat;  % from M1
load signal_rx.mat;  % from M4

%% Parameters for delay estimation:
scale = 5;
window = 100;
orders = 1;

%% Load transmitted and received signal:

[delayEst,Order,~] = MultiScale_LAP_Param(signal_tx, signal_rx, scale,window,orders);
signalEst = imshift(signal_tx,1i.*delayEst);
MSE = mean(abs((signal_rx - signalEst)).^2)/mean(signal_rx.^2);

fprintf('\t\t\t\tNo Iterations:\t MSE = %1.5f, Delay Difference = %1.4f\n', MSE, mean(abs(diff(delayEst)))/4);

delayEst(end)

% the ground truth delay ?according to the DoA provided in the metadata, I
% computed that is should be 4 samples). 
figure;
plot(signal_tx(1:100), 'r-o'); hold on;
plot(signal_rx(1:100), 'k-o'); legend("M1", "M4");






