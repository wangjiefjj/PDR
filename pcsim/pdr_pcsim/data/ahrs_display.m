clear all;
clc;

data = load('.\ahrsData.txt');

yaw = data(:, 1);
pitch = data(:, 2);
roll = data(:, 3);

figure;
plot(yaw);
titile('heading');