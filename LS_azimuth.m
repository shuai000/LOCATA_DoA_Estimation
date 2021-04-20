load('space32cm.mat');

% sampling frequency:
fs = 48000;

% distance between mics:
% d = [0.64,0.32,0.16,0.08,0.04, 0.04,0.08, 0.16, 0.32,0.64].';
d = 0.32;

% speed of sound:
c = 343;

% estimation of azimuth:


tau = 8.67/fs;
azimuth_est = asin((d\tau)*c)*180/pi


% % sampling frequency:
% fs = 48000;
%
% % distance between mics:
% % d = [0.64,0.32,0.16,0.08,0.04, 0.04,0.08, 0.16, 0.32,0.64].';
% d = [0.32,0.32,0.32,0.32,0.32,0.32].';
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
%
% figure;
% plot(azimuth_est, 'k'); hold on;
% plot(truth_angle(1, :), 'r'); hold on;
% plot(estimate_angle(1, :), 'b');