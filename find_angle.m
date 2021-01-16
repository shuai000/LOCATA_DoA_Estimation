function [az, el, r] = find_angle(h, R, p)

% p  the reference point of the array (Shuai, global, opitrack)
% h = truth.source.(src_names{src_idx}).position(:,t); % location of the source (global, opitrack)
v = -R'*(h - p); % transformed to array local coordinates (Shuai)
[az, el, r] = mycart2sph(v(1), v(2), v(3));

end