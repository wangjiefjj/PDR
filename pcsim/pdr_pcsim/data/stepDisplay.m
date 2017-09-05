clear all;
clc;

data = load('.\stepData.txt');

time = data(:, 1);
stepDet = data(:, 2);

figure;
plot(time, stepDet);
title('step det');