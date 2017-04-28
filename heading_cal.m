clear all;
clc;

data = load('.\pdr_data.log');
GNSS = 0;
SENSOR = 1;
gnss_heading_count = 0;

for i = 1:length(data)-1
    type = data(i, 1);
    time_tag = data(i, 2);
    
    if type == GNSS
        gnss_heading_count = gnss_heading_count + 1;
        if gnss_heading_count == 1
            x1_lla = [data(i,3), data(i,4), 2];
        else
            x2_lla = [data(i,3), data(i,4), 2];
            gnss_heading_location(gnss_heading_count) = heading_between_2points(x1_lla, x2_lla);
            gnss_heading_velocity(gnss_heading_count) = data(i, 9);
            gnss_heading_timing(gnss_heading_count) = time_tag;
            x1_lla = x2_lla;
        end

        if time_tag == 161781
            fprintf('time: %lu\r\n', time_tag);
            fprintf('velocity_heading: %.10f\r\n', gnss_heading_velocity(gnss_heading_count)*180/pi);
            fprintf('location_heading: %.10f\r\n', gnss_heading_location(gnss_heading_count)*180/pi);
        end
    end
end

figure;
t = gnss_heading_timing; 
y = gnss_heading_velocity*180/pi;
plot(t, y, 'b');
hold on;
y = gnss_heading_location*180/pi;
plot(t, y, 'r');