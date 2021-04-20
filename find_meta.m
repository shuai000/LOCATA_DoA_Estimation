function meta_segment = find_meta(data_segment, timestamps, time_array)

meta_segment = struct();
Ns = length(timestamps);

for k=1:length(data_segment)
    t_label = time_array(data_segment(k).index(1):data_segment(k).index(2));
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