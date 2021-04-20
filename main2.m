function [true_delay_all, evaluate_time, audio_array, audio_source, mic_positions, target, h_point, R, az_truth, el_truth] = main2(array_dir, this_array, index_all, field)

% applied for dynamic (the source is moving) case

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

evaluate_time = in_localization.timestamps;

%% Extract ground truth
%
% position_array stores all optitrack measurements.
% Extract valid measurements only (specified by required_time.valid_flag).

if is_dev
    truth = get_truth(this_array, position_array, position_source, required_time, is_dev);
else
    truth = get_truth(this_array, position_array, [], required_time, is_dev);
end

h_point = truth.array.position(:, 1);
R = squeeze(truth.array.rotation(:, 1, :));
az_truth = truth.source.(field).azimuth;
el_truth = truth.source.(field).elevation;


target = position_source.data;
%% Separate ground truth into positions of arrays (used for localization) and source position (used fo metrics)
in_localization.array = truth.array;
in_localization.array_name = this_array;
in_localization.mic_geom = truth.array.mic; % microphone position x, y, z

mic_positions = squeeze(truth.array.mic(:, 1, :));

% compute distance
source_location_all = truth.source.(field).position;
mic_location = squeeze(in_localization.mic_geom(:, 1, :));
ni = size(index_all, 1);
ns = size(source_location_all, 2);
true_delay_all = zeros(ni, ns);

for t=1:size(source_location_all, 2)
    source_location = source_location_all(:, t);
    d = sqrt(sum((source_location - mic_location) .^ 2, 1));
    delay = d ./ 343;
    delay_sample = delay .* in_localization.fs;
    
    true_delay = zeros(1, ni);
    for i=1:size(index_all, 1)
        index_tx = index_all(i, 1);
        index_rx = index_all(i, 2);
        true_delay(i) = delay_sample(index_tx) - delay_sample(index_rx);
    end
%     
%     [x, y, z] = location_search1(true_delay, mic_location, index_all);
%     [x1, y1, z1] = location_search2(true_delay, mic_location, index_all);
    
    true_delay_all(:, t) = true_delay';
end

end
