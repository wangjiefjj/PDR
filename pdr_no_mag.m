clear all;
clc;
addpath('.\functions\');

%% load SENSOR and GNSS data

% format:
% type: GNSS(0)   TimeTag(ms) Latitude(rad) Longitude(rad) Altitude(m) VelocityE(m/s2) VelocityN(m/s2) VelocityU(m/s2) Heading(rad)
% type: SENSOR(1) TimeTag(ms) AccX(m/s2) AccY(m/s2) AccZ(m/s2) GyroX(rad) GyroY(rad) GyroZ(rad) MagX(uT) MagY(uT) MagZ(uT) 
data = load('.\data\pdr_data\2017_4_28_1.log');
GNSS = 0;
SENSOR = 1;

%% main loop start
Gravity = 9.8;

% sensor heading variable
gnss_heading_time = 0;
sensor_heading_time = 0;
last_time_tag = 0;
gyro_bias = zeros(1, 3);
gyro_calibration_arry = zeros(250, 3);
acc_calibration_array = zeros(250, 1);
calibration_count = 0;

% step detection variable
step_latitude = 0;
step_longitude = 0;
ave_num = 5;
ave_count = 0;
ave_sum = 0;
acc_det_filtered = 0;
acc_array_index = 0;
acc_det_filtered_array = zeros(1,1);
acc_time_tag_array = zeros(1,1);
acc_data_process_flag = 0;
move_count = 0;
mov_window_length = 10;
acc_det_mean = 0;
pre_step_det = Gravity;
pre_slop = 0;
delta_step_time = 0;
step_lock_time = 0.3;
step_count = 0;
step_detect_flag = 0;
step_length = 0.7;

% gnss variable
step_start_flag = 0;
step_heading_flag = 0;
step_start_count = 0;
RE = 6378137.0;
esqu = 0.00669437999013;
gnss_fusion_count = 0;
gnss_heading_count = 0;
gnss_heading_window = 1;
gnss_heading_array = zeros(1,gnss_heading_window);

% step length variable
step_adjust_window = 10;
step_adjust_count = 0;

for i = 1:length(data)
    type = data(i, 1);
    time_tag = data(i, 2);
    %% gyro calibration (5s static data)
    if type == SENSOR
        acc_data = data(i, 3:5);
        gyro_data = data(i, 6:8);
        acc_det = sqrt(sum(acc_data.^2));
        calibration_count = calibration_count + 1;
        if calibration_count > 250
            for j = 2 : 250
                gyro_calibration_arry(j-1, :) = gyro_calibration_arry(j, :);
                gyro_calibration_arry(250, :) = gyro_data;
                acc_calibration_array(j-1, :) = acc_calibration_array(j, :);
                acc_calibration_array(250, :) = acc_det;
            end
        else
            gyro_calibration_arry(calibration_count, :) = gyro_data;
            acc_calibration_array(calibration_count) = acc_det;
        end

        if calibration_count == 250
            acc_energy = mean(acc_calibration_array);
            acc_std = std(acc_calibration_array);
            if acc_std < 0.1
                gyro_bias = mean(gyro_calibration_arry)
                fprintf('static detected and start gyro calibration at time %lu\r\n', time_tag);
            end
        end
    end
        
    if type == SENSOR && step_start_flag ~= 0   
        acc_data = data(i, 3:5);
        gyro_data = data(i, 6:8);
        acc_det = sqrt(sum(acc_data.^2));
        
       %% smooth acc data
        ave_count = ave_count + 1;
        if ave_count == ave_num
            ave_count = 0;
            acc_array_index = acc_array_index + 1;
            acc_det_filtered = (ave_sum + acc_det)/ave_num;
            ave_sum = 0;
            acc_det_filtered_array(acc_array_index) = acc_det_filtered;
            acc_time_tag_array(acc_array_index) = time_tag/1000;
            acc_data_process_flag = 1;
            step_det = acc_det_filtered;
        else
            ave_sum = ave_sum + acc_det;
        end
        
       %% heading calculate (50Hz)
        if gnss_heading_time ~= 0
            % caluculate heading by gyro
            if abs(gyro_data(3))*180/pi < 0.01
%                 gyro_data(1) = 0;
%                 gyro_data(2) = 0;
%                 gyro_data(3) = 0;
%                 fprintf('gyro heading not compute at time: %lu\r\n', time_tag);
            else
%                 fprintf('gyro heading compute at time: %lu\r\n', time_tag);
            end
            heading_aid_inteval = (time_tag - gnss_heading_time)/1000;   % s
            if heading_aid_inteval < 30
                dt = (time_tag - sensor_heading_time)/1000;
                sensor_heading_time = time_tag;
                q = zeros(4, 1);
                q = euler2q(step_heading, 0, 0);
                Cbn = q2dcm(q);
                Cnb = Cbn';
                Wepp = zeros(3, 1); % latitude and longitude rate
                Wiep = zeros(3, 1); % earth rate in the navigation frame
                Wipp = Wiep + Wepp;
                Wipb = Cnb * Wipp;
                Wpbb = gyro_data' - gyro_bias' - Wipb;

                dq = zeros(4, 1);
                dq(1) = -(Wpbb(1)*q(2) + Wpbb(2)*q(3) + Wpbb(3)*q(4))/2;
                dq(2) = (Wpbb(1)*q(1) + Wpbb(3)*q(3) - Wpbb(2)*q(4))/2;
                dq(3) = (Wpbb(2)*q(1) - Wpbb(3)*q(2) + Wpbb(1)*q(4))/2;
                dq(4) = (Wpbb(3)*q(1) + Wpbb(2)*q(2) - Wpbb(1)*q(3))/2;

                q = q + dq*dt;
                q = q_norm(q);
                Cbn = q2dcm(q);
                [yaw, pitch, roll] = dcm2euler(Cbn);
                step_heading = yaw;
            else
                fprintf('PDR has no heading aided at time: %lu, delta time: %.10fs\r\n', time_tag, heading_aid_inteval);
            end
        end
        
        %% step detecion (10Hz)
         if acc_data_process_flag
             acc_data_process_flag = 0;
             % det move average calculate for threshold estimate and the window length is 10 (1s)
             if move_count < mov_window_length
                 move_count = move_count + 1;
                 acc_det_mean = (acc_det_mean*(move_count-1) + step_det)/move_count;
             else 
                 acc_det_mean = (acc_det_mean*(mov_window_length-1) + step_det)/mov_window_length;
             end
             
              % determin the threshold
             if abs(acc_det_mean - Gravity) > 0.8
                 threshold = 0.15 * acc_det_mean;
             else
                 acc_det_mean = Gravity;
                 threshold = 0.05 * acc_det_mean;
             end
             
             % slop calculate
            if step_det - pre_step_det > 0
                slop = 1;
            else
                slop = -1;
            end
            
            % delta step time calculate
            if acc_array_index == 1
                delta_step_time = 0.1;
            else
                delta_step_time = delta_step_time + acc_time_tag_array(acc_array_index) - acc_time_tag_array(acc_array_index-1);
            end
            
            % step detect (wave summit & delta time & step det)
            if ( pre_slop > 0 && slop < 0  && delta_step_time >= step_lock_time && pre_step_det - acc_det_mean > threshold)
                step_count = step_count + 1;
                step_adjust_count = step_adjust_count + 1;
                step_detect_flag = 1;
                step_timing(step_count) = acc_time_tag_array(acc_array_index-1);
                delta_step_time = 0;
            else
                if delta_step_time > 3 % 3s has no step for a long time, need to reset. 
                   delta_step_time = 0;
                end
            end

            pre_slop = slop;
            pre_step_det = step_det; 
            
           %% PDR estimate latitude and longitude 
            if step_detect_flag ~= 0
                step_detect_flag = 0;
                
                % calculate the step trajectory
                RM = RE * (1 - esqu) / ((1 - esqu * sin(step_latitude) * sin(step_latitude))^1.5);
                RN = RE / (sqrt(1 - esqu * sin(step_latitude) * sin(step_latitude)));

                step_latitude = step_latitude + step_length*cos(step_heading)/RM;
                step_longitude = step_longitude + step_length*sin(step_heading)/(RN*cos(step_latitude));
                step_altitude = gnss_altitude;
                
                % store all step trajectory
                step_latitude_array(step_count) = step_latitude*180/pi;
                step_longitude_array(step_count) = step_longitude*180/pi;
                step_altitude_array(step_count) = step_altitude;
            end
             
         end
    elseif type == GNSS
       %% gnss data process
        gnss_latitude = data(i, 3);
        gnss_longitude = data(i, 4);
        gnss_altitude = data(i, 5);
        gnss_vel_ENU = data(i, 6:8);
        % determin the start timing for PDR
        % step velocity > 0.5m/s and keep 10s means step start.
        gnss_vel = sqrt(sum(gnss_vel_ENU.^2));
        if step_start_flag == 0
            if gnss_vel > 0.5
                step_start_count = step_start_count + 1;
                if step_start_count >= 10
                    step_start_flag = 1;
                    step_latitude = gnss_latitude;
                    step_longitude = gnss_longitude;
                    step_heading = data(i, 9);
                    fprintf('PDR start at time: %lu\r\n', time_tag);
                end
            else
                step_start_count = 0;
            end
        end
        % check gnss heading in 5 seconds window
%         if step_start_flag
%             if gnss_vel > 0.5
%                 gnss_heading = data(i, 9);
%                 gnss_heading_count = gnss_heading_count + 1;
%                 if gnss_heading_count <= gnss_heading_window
%                     gnss_heading_array(gnss_heading_count) = gnss_heading;
%                 else
%                     for j = 2:gnss_heading_window
%                         gnss_heading_array(j-1) = gnss_heading_array(j);
%                     end
%                     gnss_heading_array(gnss_heading_window) = gnss_heading;
%                     gnss_heading_mean = mean(gnss_heading_array);
%                     gnss_heading_std = std(gnss_heading_array);
%                     if gnss_heading_std < 5/180*pi % 5 deg
%                         step_heading = gnss_heading_mean;
%                         gnss_heading_time = time_tag;
% %                         fprintf('Velocity Heading availiale at time:%lu  Heading value is:%.10f\r\n', time_tag, gnss_heading_mean);
%                     end
%                 end
%             else
%                 gnss_heading_count = 0;
%             end
%         end
        
        % low speed condition, the gnss velocity and heading is not accurate, and use location estimating the heading
%         if step_start_flag
%             if gnss_vel > 0
%                 if step_heading_flag == 0
%                     x1_lla = [data(i,3), data(i,4), 2];
%                     gnss_heading = data(i, 9);
%                     step_heading_flag = 1;
%                 else
%                     x2_lla = [data(i,3), data(i,4), 2];
%                     gnss_heading = heading_between_2points(x1_lla, x2_lla);
%                     x1_lla = x2_lla;
%                 end
%                 step_heading = gnss_heading;
%                 gnss_heading_time = time_tag;
%                 sensor_heading_time = time_tag;
%             end
%         end

        if step_start_flag
            if step_heading_flag == 0
                x1_lla = [data(i,3), data(i,4), 2];
                gnss_heading = data(i, 9);
                step_heading_flag = 1;
            else
                x2_lla = [data(i,3), data(i,4), 2];
                gnss_heading = heading_between_2points(x1_lla, x2_lla);
                x1_lla = x2_lla;
            end
            
            % south heading discontinued
            if gnss_heading > -pi && gnss_heading < -pi + 2/180*pi
                fprintf('South heading at time:%lu  Heading value is:%.10f\r\n', time_tag, gnss_heading);
                gnss_heading = gnss_heading + 2*pi;
            end
            
            if gnss_heading > pi - 2/180*pi && gnss_heading < pi
                fprintf('South heading at time:%lu  Heading value is:%.10f\r\n', time_tag, gnss_heading);
                gnss_heading = gnss_heading;
            end
            
            % location heading filter
            if step_heading_flag == 1
                if gnss_vel > 0
                    gnss_heading_count = gnss_heading_count + 1;
                    if gnss_heading_count <= gnss_heading_window
                        gnss_heading_array(gnss_heading_count) = gnss_heading;
                    else
                        for j = 2:gnss_heading_window
                            gnss_heading_array(j-1) = gnss_heading_array(j);
                        end
                        gnss_heading_array(gnss_heading_window) = gnss_heading;
                        gnss_heading_mean = mean(gnss_heading_array);
                        gnss_heading_std = std(gnss_heading_array);
                        if gnss_heading_std < 20/180*pi % 20 deg
                            step_heading = gnss_heading_mean;
                            gnss_heading_time = time_tag;
                            sensor_heading_time = time_tag;
%                             fprintf('Heading availiale at time:%lu  Heading value is:%.10f\r\n', time_tag, gnss_heading_mean);
                        end
                    end
                else
                    gnss_heading_count = 0;
                end
            end
        end
        
%         if gnss_vel > 1.0
%             step_heading = data(i, 9);
%             gnss_heading_time = time_tag;
%             sensor_heading_time = time_tag;
%         end
        
       %% fusion pdr trajectory and gnss position every 10 seconds
        gnss_fusion_count = gnss_fusion_count + 1;
        if gnss_fusion_count >= 10
            gnss_fusion_count = 0;
            step_latitude = 0.35*gnss_latitude + 0.65*step_latitude;
            step_longitude = 0.35*gnss_longitude + 0.65*step_longitude;
        end
        
       %% adaptive step length process (every 20 step, about 10s)
        if step_adjust_count > step_adjust_window
            delta_latitude = gnss_latitude - adjust_latitude_start;
            delta_longitude = gnss_longitude - adjust_longitude_start;
            ave_latitude = (gnss_latitude + adjust_latitude_start)/2;
            delta_N = delta_latitude * RM;
            delta_E = delta_longitude * RN * cos(ave_latitude);
            distance = sqrt(delta_N.^2 + delta_E.^2);
            adjust_length = distance / step_adjust_count;
            step_adjust_count = 0;
            if adjust_length > 0.1 && adjust_length < 2.0
                step_length = adjust_length;
            end
        end
        if step_adjust_count == 0
            adjust_latitude_start = gnss_latitude;
            adjust_longitude_start = gnss_longitude;
        end

    end
    
end

%% plot the filter acc det data for debug
Fs = 10;                            % About 10Hz
L = length(acc_det_filtered_array); % Length of signal
t = acc_time_tag_array;             % Time vector
y = acc_det_filtered_array;

if 0
figure;
start_time = 600*Fs + 1;
sample_time = 100;                  % Display "sample_time" length data
plot(t(start_time:start_time + Fs*sample_time),y(start_time:start_time + Fs*sample_time))
title('Filtered Acc Det in 10 Second')
xlabel('Time (seconds)')
ylabel('Y(t)')
end

%% plot the filter gyro data for debug
if 0
figure;
start_time = 0*Fs + 1;
sample_time = 100;                  % Display "sample_time" length data
y = gyro_filtered_array(:,1)*180/pi;       % x
plot(t(start_time:start_time + Fs*sample_time),y(start_time:start_time + Fs*sample_time), 'r')
hold on;
y = gyro_filtered_array(:,2)*180/pi;       % y
plot(t(start_time:start_time + Fs*sample_time),y(start_time:start_time + Fs*sample_time), 'g')
y = gyro_filtered_array(:,3)*180/pi;       % z
plot(t(start_time:start_time + Fs*sample_time),y(start_time:start_time + Fs*sample_time), 'b')
title('Filtered Gyro Data')
xlabel('Time (seconds)')
ylabel('Value (Deg/s)')
legend('x axis', 'y axis', 'z axis')
end

%% convert pdr trajectory to KML file
gps2kml('gpsOut',step_timing,step_latitude_array,step_longitude_array,step_altitude_array,'o-b','MarkerSize',10,'LineWidth',3);
