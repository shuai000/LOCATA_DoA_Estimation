% Test angle

clc; clear all; close all; 

main('D:\LOCATA\dev\', 'D:\LOCATA\results\', 1, {'dicit'}, [3]);

x = -1;
y = -1;
z = -1;

[az, el, r] = mycart2sph(x, y, z);

rad2deg(az)
rad2deg(el)


