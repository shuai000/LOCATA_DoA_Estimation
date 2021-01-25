function [data_struct, meta_struct] = find_segment(data, labels, t_array, timestamps)
%%% Return the signal for each segment based on the active labels
%%% data_struct: struct array, field: data, index

data_struct = struct();

indicator = diff(labels);

index_start = find(indicator == 1);
index_end = find(indicator == -1);

N = min(length(index_start), length(index_end));

for i=1:N
    s_i = index_start(i);
    e_i = index_end(i);
    
    data_struct(i).index_interval = [s_i, e_i];
    data_struct(i).data = data(s_i:e_i, :);
    data_struct(i).time = t_array(s_i:e_i);
end

if length(index_start) - length(index_end) == 1
    s_i = index_start(end);
    e_i = size(data, 1);
    
    data_struct(i+1).index_interval = [s_i, e_i];
    data_struct(i+1).data = data(s_i:e_i, :);
    data_struct(i+1).time = t_array(s_i:e_i);
end

meta_struct = find_meta(data_struct, timestamps, t_array);

end

function meta_segment = find_meta(data_segment, timestamps, time_array)

meta_segment = struct();
Ns = length(timestamps);

for k=1:length(data_segment)
    t_label = time_array(data_segment(k).index_interval(1):data_segment(k).index_interval(2));
    index_array = [];
    true_array = [];
    for i = 1:Ns
        % find index to extract delay estimation, not perfertly synchronized
        current_index = find_index(timestamps(i), t_label);
        if current_index > 0
            index_array = [index_array, current_index]; % w.r.t. global data index, (long)
            true_array = [true_array, i]; % w.r.t. timestamp (short)
        end
    end
    meta_segment(k).local_index = index_array;
    meta_segment(k).true_index = true_array;
    
end

end