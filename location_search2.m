
function [x, y, z, min_error] = location_search2(estimated_delay, mic_positions, index_array)
% IN GLOBAL COORDINATE, DEFINED BY THE OPiTRACK SYSTEM

x = 0;
y = 0;
z = 0;

my_step = .1;
min_error = 1000000;
error_array = [];

for xc = -3:my_step:3
    for yc = -.1:my_step:4
        for zc = 1.3:my_step:1.75
            source_cartisan = [xc, yc, zc]';
            % compute distance
            dis_ms = sqrt(sum((source_cartisan - mic_positions) .^ 2, 1));
            
            % compute delay
            v_speed = 343;
            dt = dis_ms ./ v_speed * 48000;  % by samples
            delay_anticipate = zeros(1, size(index_array, 1));
            for i = 1:size(index_array, 1)
                delay_anticipate(i) = dt(index_array(i, 1)) - dt(index_array(i, 2));
            end
            error_delay = sqrt(sum((estimated_delay - delay_anticipate) .^ 2));
            
            if error_delay < min_error
                x = xc; y = yc; z = zc;
                min_error = error_delay;
            end
            
            error_array = [error_array, error_delay];
            
        end
    end
end

% figure;
% plot(error_array);
end