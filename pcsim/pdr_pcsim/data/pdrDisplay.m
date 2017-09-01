clear all;
clc;

data = load('.\output.txt');

time = data(:, 1);
latitude = data(:, 2);
longitude = data(:, 3);
altitude = data(:, 4);

gps2kml('pdrOut', time, latitude, longitude, altitude, 'o-b', 'MarkerSize', 10, 'LineWidth', 3);