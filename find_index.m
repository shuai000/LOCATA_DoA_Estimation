function index = find_index(current_time, t_array)
index = 0;
N = length(t_array);

for i=1:N-1
    if current_time >= t_array(i) && current_time < t_array(i+1)
        index = i;
        break;
    end
end