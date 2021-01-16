function true_delay = main1(array_dir, this_array, index_all)

% function main_eval(data_dir, results_dir, arrays, tasks)
% main function for participants of the LOCATA Challenge to load each recording
% for specific tasks and arrays to run their algorithm for the Eval or Dev database.
%
% Inputs:
%   data_dir:   String with directory path for the LOCATA Dev or Eval database
%   save_dir:   String with directory path in which to save the results of this
%               function
%   is_dev:     Kind of database specified by data_dir
%               0: Eval database
%               1: Dev database
%   arrays:     Cell array with array names which should be evaluated (optional)
%               Cell array {'benchmark2', 'eigenmike', 'dicit','dummy'} is taken
%               as default which contains all available arrays
%   tasks:      Vector with task(s) (optional)
%               Vector [1,2,3,4,5,6] is taken as default which evaluates
%               over all available tasks
%
% Outputs: N/A (saves results as csv files in save_dir)
%
% Authors: Christine Evers, c.evers@imperial.ac.uk
%          Heiner Loellmann, Loellmann@LNT.de
%
% Reference: LOCATA documentation for participants (v3)
%            www.locata-challenge.org
%
% Notice: This programm is part of the LOCATA evaluation release.
%         Please report problems and bugs to info@locata-challenge.org.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add matlab directory and sub-folders to path:
addpath(genpath('./'))


%% Load data

is_dev = 1;
if is_dev
    [audio_array, audio_source, position_array, position_source, required_time] = load_data(array_dir, is_dev);
else
    [audio_array, position_array, required_time] = load_data(array_dir, is_dev);
end

%% Load signal

% Get number of mics and mic array geometry:
in_localization.numMics = size(position_array.data.(this_array).mic,3);

% Signal and sampling frequency:
in_localization.y = audio_array.data.(this_array)';      % signal
in_localization.fs = audio_array.fs;                     % sampling freq

duration_time = size(in_localization.y, 2) / in_localization.fs;

%% Users must provide estimates for each time stamp in in.timestamps

% Time stamps required for evaluation
in_localization.timestamps = elapsed_time(required_time.time);
in_localization.timestamps = in_localization.timestamps(find(required_time.valid_flag));
in_localization.time = required_time.time(:,find(required_time.valid_flag));

%% Extract ground truth
%
% position_array stores all optitrack measurements.
% Extract valid measurements only (specified by required_time.valid_flag).

if is_dev
    truth = get_truth(this_array, position_array, position_source, required_time, is_dev);
else
    truth = get_truth(this_array, position_array, [], required_time, is_dev);
end

%% Separate ground truth into positions of arrays (used for localization) and source position (used fo metrics)

in_localization.array = truth.array;
in_localization.array_name = this_array;
in_localization.mic_geom = truth.array.mic; % microphone position x, y, z

% compute distance
source_location = truth.source.loudspeaker3.position(:, 1); % static, take the first
mic_location = squeeze(in_localization.mic_geom(:, 1, :));

d = sqrt(sum((source_location - mic_location) .^ 2, 1));
delay = d ./ 343;
delay_sample = delay .* in_localization.fs;

true_delay = zeros(1, size(index_all, 1));
for i=1:size(index_all, 1)
    index_tx = index_all(i, 1);
    index_rx = index_all(i, 2);
    true_delay(i) = delay_sample(index_tx) - delay_sample(index_rx);
end

%             % test
%             a1 = 1;
%             a2 = 4;
%
%             p1 = mic_location(:, a1);
%             p2 = mic_location(:, a2);
%             sqrt(sum((p1 - p2) .^2))


end
