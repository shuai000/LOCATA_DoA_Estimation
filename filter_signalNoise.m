function [y, index_y] = filter_signalNoise(x, threshold)
N = length(x);
y = [];
index_y = [];

for i=1:N
    if abs(x(i)) > threshold
        y = [y, x(i)];
        index_y = [index_y, i];
    end
end

end