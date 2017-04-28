clear all;
clc;
addpath('.\functions\');

data = load('.\pdr_data.log');
GNSS = 0;
SENSOR = 1;
count = 0;

for i = 1:length(data)
    type = data(i, 1);
    time_tag = data(i, 2);
    if type == SENSOR
        count = count + 1;
        gyro(count, :) = data(i, 6:8);
        time(count, :) = time_tag/1000;
    end
end

plot(time, gyro(:, 1), 'r');
hold on;
plot(time, gyro(:, 2), 'g');
plot(time, gyro(:, 3), 'b');

x_offset = mean(gyro(1:250, 1))
y_offset = mean(gyro(1:250, 2))
z_offset = mean(gyro(1:250, 3))