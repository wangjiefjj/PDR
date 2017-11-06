clear all;
clc;
addpath('.\functions\basic');
addpath('.\functions\mag_calibration');
addpath('.\functions\misc');

%% load SENSOR and GNSS data

% format:
% type: GNSS(0)   TimeTag(ms) Latitude(rad) Longitude(rad) Altitude(m) VelocityE(m/s2) VelocityN(m/s2) VelocityU(m/s2) Heading(rad)
% type: SENSOR(1) TimeTag(ms) AccX(m/s2) AccY(m/s2) AccZ(m/s2) GyroX(rad) GyroY(rad) GyroZ(rad) MagX(uT) MagY(uT) MagZ(uT) 
data = load('.\data\pdr_data.log');
GNSS = 0;
SENSOR = 1;

%% configure

SampleRate = 50;   % Hz
AlignmentTime = 2; % second
MagCalibrationTime = 10;   % second


%% prepare data for mag calibration & initial alignment
count = 1;
for i = 1:length(data)
    type = data(i, 1);
    if type == SENSOR
        Mag_Raw(count, :) = data(i, 9:11);
        Acc(count, :) = data(i, 3:5);
        Gyro(count, :) = data(i, 6:8);
        count = count + 1;
    end
end

N = SampleRate*MagCalibrationTime;
[Cali_B, Cali_V, Cali_W_inv, Cali_Error] = cali7eig(Mag_Raw(1:N, :));
if Cali_Error < 0.05
    mag_calibration = 1;
else
    mag_calibration = 0;
    disp('mag calibration can not be execuated');
end

%% initial alignment

% find the first position fix index
count = 1;
for i = 1:length(data)
    type = data(i, 1);
    if type == SENSOR
        count = count + 1;
    else
        break;
    end
end
M = count;
N = SampleRate*AlignmentTime;
initial_alignment_flag = 0;
acc_alignment = Acc(M:M+N, :);
mag_alignment = Mag_Raw(M:M+N, :);

% calibrated mag data
for i = 1 : length(mag_alignment)
    Bp = mag_alignment(i, :);
    Bc = Cali_W_inv*(Bp - Cali_V)';
    mag_alignment(i, :) = Bc';
end

% static check
acc_x_var = var(acc_alignment(1:N, 1));
acc_y_var = var(acc_alignment(1:N, 2));
acc_z_var = var(acc_alignment(1:N, 3));
mag_x_var = var(mag_alignment(1:N, 1));
mag_y_var = var(mag_alignment(1:N, 2));
mag_z_var = var(mag_alignment(1:N, 3));
acc_var_threshold = 0.01;
mag_var_threshold = 1;
if acc_x_var < acc_var_threshold && acc_y_var < acc_var_threshold && acc_z_var < acc_var_threshold && ...
   mag_x_var < mag_var_threshold && mag_y_var < mag_var_threshold && mag_z_var < mag_var_threshold
   initial_alignment_flag = 1;
end

if initial_alignment_flag == 1
    g_x_mean = -mean(acc_alignment(1:N, 1));
    g_y_mean = -mean(acc_alignment(1:N, 2));
    g_z_mean = -mean(acc_alignment(1:N, 3));
    mag_x_mean = mean(mag_alignment(1:N, 1));
    mag_y_mean = mean(mag_alignment(1:N, 2));
    mag_z_mean = mean(mag_alignment(1:N, 3));
    Cnb_initial = ecompass_ned([g_x_mean, g_y_mean, g_z_mean], [mag_x_mean, mag_y_mean, mag_z_mean]);
    Cbn_initial = Cnb_initial';
    [yaw_initial, pitch_initial, roll_initial] = dcm2euler(Cbn_initial);
    q = euler2q(yaw_initial, pitch_initial, roll_initial);
    q_pure = q;
    % compute the geomagnetic inclination angle
    geoB = norm([mag_x_mean, mag_y_mean, mag_z_mean]);
    gmod = norm([g_x_mean, g_y_mean, g_z_mean]);
    geo_inclination = asin(dot([mag_x_mean, mag_y_mean, mag_z_mean], [g_x_mean, g_y_mean, g_z_mean])/geoB/gmod); % rad
    if mag_calibration == 1
        Mag_vector = [Cali_B*cos(geo_inclination), 0, Cali_B*sin(geo_inclination)]';
    else
        Mag_vector = [geoB*cos(geo_inclination), 0, geoB*sin(geo_inclination)]';
        Cali_B = geoB;
    end
end

%% general variable
Ge = 9.80665;
RE = 6378137.0;
esqu = 0.00669437999013;
G_vector = [0, 0, Ge]';
sensor_count = 0;
static_count = 0;
initial_latitude = 0;
initial_longitude = 0;
initial_altitude = 0;

%% ahrs variable
gyro_bias = zeros(3, 1);
acc_bias = zeros(3, 1);
gyro_smooth_count = 0;
gyro_smooth_length = 50;
gyro_smooth_array = zeros(gyro_smooth_length, 3);

x = zeros(9, 1); % roll, pitch, yaw, gyro_bias_x, gyro_bias_y, gyro_bias_z, acc_bias_x, acc_bias_y, acc_bias_z
F = zeros(9, 9);
PHIM = zeros(9, 9);
qdt = zeros(9, 9);
Q = zeros(9, 9);
G = zeros(9, 9);
Corr_time_gyro = 10;
Corr_time_acc = 10;
sigma_Win = 1.0e-6;
sigma_acc = ((5.0e-4) * 9.78032667 * (5.0e-4) * 9.78032667);
sigma_gyro = (2 * pi / 180.0 / 3600 * 2 * pi / 180.0 / 3600);
sigma_phim_e_n = 1.0*pi/180;
sigma_phim_u = 1.0*pi/180;
sigma_phim_gyro = 10*pi/180/3600;
sigma_phim_acc = 0.3;
P = eye(9, 9);
P(1, 1) = sigma_phim_e_n^2;
P(2, 2) = sigma_phim_e_n^2;
P(3, 3) = sigma_phim_u^2;
P(4, 4) = sigma_phim_gyro^2;
P(5, 5) = sigma_phim_gyro^2;
P(6, 6) = sigma_phim_gyro^2;
P(7, 7) = sigma_phim_acc^2;
P(8, 8) = sigma_phim_acc^2;
P(9, 9) = sigma_phim_acc^2;

%% step detection variable
step_latitude = 0;
step_longitude = 0;
ave_num = 5;
ave_count = 0;
ave_sum = 0;
acc_det_filtered = 0;
acc_array_index = 0;
acc_data_process_flag = 0;
move_count = 0;
mov_window_length = 10;
acc_det_mean = 0;
pre_step_det = Ge;
pre_slop = 0;
delta_step_time = 0;
step_lock_time = 0.3;
step_count = 0;
step_detect_flag = 0;
step_length = 0.8;
step_freq = 0;
step_start_flag = 0;

%% pdr gnss fusion variable
Corr_time_step_length = 100;
sigma_n = 0.1;
sigma_e = 0.1;
sigma_step_length = 0.1;
sigma_step_heading = 5/180*pi;
sigma_p = [0.1, 0.1, 0.1, 5/180*pi];
pdr_p = diag(sigma_p.^2);
pdr_x = zeros(4, 1);

%% main loop
for i = M:length(data)
    type = data(i, 1);
    time_tag = data(i, 2);
    %% GNSS aiding
    if type == GNSS
        gnss_latitude = data(i, 3);
        gnss_longitude = data(i, 4);
        gnss_altitude = data(i, 5);
        gnss_vel = data(i, 6:8);
        gnss_heading = data(i, 9);
        if step_start_flag ~= 1
            static_count = static_count + 1;
            initial_latitude = initial_latitude + gnss_latitude;
            initial_longitude = initial_longitude + gnss_longitude;
            initial_altitude = initial_altitude + gnss_altitude;
            
            if static_count == 10
                step_start_flag = 1;
                pdr_latitude = initial_latitude / 10;
                pdr_longitude = initial_longitude / 10;
                pdr_altitude = initial_altitude / 10;
            else
                continue;
            end
        end
        
        gnss_vel_det = norm(gnss_vel);
        if gnss_vel_det > 1.5
            dt = 1;
            
            pdr_f = zeros(4, 4);
            pdr_f(1, 3) = step_freq*cos(yaw);
            pdr_f(1, 4) = -step_length*step_freq*sin(yaw);
            pdr_f(2, 3) = step_freq*sin(yaw);
            pdr_f(2, 4) = step_length*step_freq*cos(yaw);
            pdr_f(3, 3) = -1/Corr_time_step_length;
            pdr_f(4, 4) = 0;
            
            pdr_qdt = zeros(4, 4);
            pdr_qdt(1, 1) = sigma_n^2;
            pdr_qdt(2, 2) = sigma_e^2;
            pdr_qdt(3, 3) = sigma_step_length^2;
            pdr_qdt(4, 4) = sigma_step_heading^2;
            
            pdr_g = eye(4, 4);
            
            % Q matrix discretization-2 order
            Q_basic = pdr_g*pdr_qdt*pdr_g';
            M1 = Q_basic;
            M2 = Q_basic*pdr_f'+pdr_f*Q_basic;
            pdr_q = dt*M1 + 1/2*dt*dt*M2;

            % PHIM matrix discretization-2 order
            I = eye(4, 4);
            pdr_phim = I + dt*pdr_f + 1/2*dt*dt*pdr_f*pdr_f;        

            % predict
            pdr_x = pdr_phim*pdr_x;
            pdr_p = pdr_phim*pdr_p*pdr_phim' + pdr_q;

            % update from gnss
            gnss_N = gnss_latitude * RM;
            gnss_E = gnss_longitude * RN;
            pdr_N = pdr_latitude * RM;
            pdr_E = pdr_longitude * RN;
            if gnss_heading - yaw > pi
                gnss_heading = gnss_heading - 2*pi;
            end
            if gnss_heading - yaw < -pi
                gnss_headng = gnss_heading + 2*pi;
            end
            pdr_z = [gnss_N - pdr_N; gnss_E - pdr_E; gnss_heading - yaw];
            pdr_h = zeros(3, 4);
            pdr_h(1, 1) = 1;
            pdr_h(2, 2) = 1;
            pdr_h(3, 4) = 1;
            sigma_r = [20, 20, 10/180*pi];
            pdr_r = diag(sigma_r.^2);
            pdr_i = eye(4, 4);
            pdr_k = pdr_p*pdr_h'*((pdr_h*pdr_p*pdr_h'+pdr_r)^-1);
            pdr_x = pdr_x + pdr_k*(pdr_z - pdr_h*pdr_x);
            pdr_p = (pdr_i - pdr_k*pdr_h)*pdr_p;
            
            if pdr_x(1) < 0.5
                pdr_N = pdr_N + pdr_x(1);
                pdr_latitude = pdr_N / RM;
            end
            pdr_x(1) = 0;
            
            if pdr_x(2) < 0.5
                pdr_E = pdr_E + pdr_x(2);
                pdr_longitude = pdr_E / RN;
            end
            pdr_x(2) = 0;

            if abs(pdr_x(4)) < 10*pi/180
                yaw = yaw + pdr_x(4);
                if yaw > pi
                    yaw = yaw - 2*pi;
                end
                if yaw < -pi
                    yaw = yaw + 2*pi;
                end
                pdr_x(4) = 0;
                q = euler2q(yaw, pitch, roll);
            end
            
            step_length = step_length + pdr_x(3);
            pdr_x(3) = 0;
            
            pdr_latitude_array(step_count) = pdr_latitude*180/pi;
            pdr_longitude_array(step_count) = pdr_longitude*180/pi;
            pdr_altitude_array(step_count) = pdr_altitude;
        end
    
    % SENSOR
    elseif type == SENSOR
        sensor_count = sensor_count + 1;
        heading_gnss(sensor_count) = gnss_heading;
        
        %% AHRS
        Mag = data(i, 9:11);
        Acc = data(i, 3:5);
        Gyro = data(i, 6:8);
        % smooth gyro data
        gyro_smooth_count = gyro_smooth_count + 1;
        if gyro_smooth_count > gyro_smooth_length
            for j = 1 : gyro_smooth_length - 1
                gyro_smooth_count = gyro_smooth_length;
                gyro_smooth_array(j, :) = gyro_smooth_array(j+1, :);
            end
            gyro_smooth_array(gyro_smooth_length, :) = Gyro; 
        else
            gyro_smooth_array(gyro_smooth_count, :) = Gyro;
        end
        slot_number = min(gyro_smooth_count, gyro_smooth_length);
        Gyro(1) = mean(gyro_smooth_array(1:slot_number, 1));
        Gyro(2) = mean(gyro_smooth_array(1:slot_number, 2));
        Gyro(3) = mean(gyro_smooth_array(1:slot_number, 3));
        for j = 1:3
            if abs(Gyro(j)) < 0.08
               Gyro(j) = 0; 
            end
        end
        
        %% pure gyro estimate
        
        Cbn = q2dcm(q_pure);
        Cnb = Cbn';

        Wpbb = Gyro';

        dq = zeros(4, 1);
        dq(1) = -(Wpbb(1)*q_pure(2) + Wpbb(2)*q_pure(3) + Wpbb(3)*q_pure(4))/2;
        dq(2) = (Wpbb(1)*q_pure(1) + Wpbb(3)*q_pure(3) - Wpbb(2)*q_pure(4))/2;
        dq(3) = (Wpbb(2)*q_pure(1) - Wpbb(3)*q_pure(2) + Wpbb(1)*q_pure(4))/2;
        dq(4) = (Wpbb(3)*q_pure(1) + Wpbb(2)*q_pure(2) - Wpbb(1)*q_pure(3))/2;
        dt = 1/SampleRate;
        q_pure = q_pure + dq*dt;
        q_pure = q_norm(q_pure);
        Cbn = q2dcm(q_pure);
        [yaw_pure, pitch_pure, roll_pure] = dcm2euler(Cbn); % pure gyro estimation
        heading_pure(sensor_count) = yaw_pure;
        
        %% fusion estimate

        % quaternion integration
        Cbn = q2dcm(q);
        Cnb = Cbn';

        Wepp = zeros(3, 1); % no latitude information in computing latitude and longitude rate
        Wiep = zeros(3, 1); % no latitude information in computing earth rate in the navigation frame
        Wipp = Wiep + Wepp;
        Wipb = Cnb * Wipp;
        Wpbb = Gyro' - gyro_bias - Wipb;

        dq = zeros(4, 1);
        dq(1) = -(Wpbb(1)*q(2) + Wpbb(2)*q(3) + Wpbb(3)*q(4))/2;
        dq(2) = (Wpbb(1)*q(1) + Wpbb(3)*q(3) - Wpbb(2)*q(4))/2;
        dq(3) = (Wpbb(2)*q(1) - Wpbb(3)*q(2) + Wpbb(1)*q(4))/2;
        dq(4) = (Wpbb(3)*q(1) + Wpbb(2)*q(2) - Wpbb(1)*q(3))/2;
        dt = 1/SampleRate;
        q = q + dq*dt;
        q = q_norm(q);
        Cbn = q2dcm(q);
        
        %% kalman filter
        F(1, 4) = -Cbn(1, 1);
        F(1, 5) = -Cbn(1, 2);
        F(1, 6) = -Cbn(1, 3);
        F(2, 4) = -Cbn(2, 1);
        F(2, 5) = -Cbn(2, 2);
        F(2, 6) = -Cbn(2, 3);
        F(3, 4) = -Cbn(3, 1);
        F(3, 5) = -Cbn(3, 2);
        F(3, 6) = -Cbn(3, 3);
        F(4, 4) = -1/Corr_time_gyro;
        F(5, 5) = -1/Corr_time_gyro;
        F(6, 6) = -1/Corr_time_gyro;
        F(7, 7) = -1/Corr_time_acc;
        F(8, 8) = -1/Corr_time_acc;
        F(9, 9) = -1/Corr_time_acc;

        qdt(1, 1) = sigma_Win;
        qdt(2, 2) = sigma_Win;
        qdt(3, 3) = sigma_Win;
        qdt(4, 4) = sigma_gyro;
        qdt(5, 5) = sigma_gyro;
        qdt(6, 6) = sigma_gyro;
        qdt(7, 7) = sigma_acc;
        qdt(8, 8) = sigma_acc;
        qdt(9, 9) = sigma_acc;

        G(1, 1) = -Cbn(1, 1);
        G(1, 2) = -Cbn(1, 2);
        G(1, 3) = -Cbn(1, 3);
        G(2, 1) = -Cbn(2, 1);
        G(2, 2) = -Cbn(2, 2);
        G(2, 3) = -Cbn(2, 3);
        G(3, 1) = -Cbn(3, 1);
        G(3, 2) = -Cbn(3, 2);
        G(3, 3) = -Cbn(3, 3);
        G(4, 4) = 1;
        G(5, 5) = 1;
        G(6, 6) = 1;
        G(7, 7) = 1;
        G(8, 8) = 1;
        G(9, 9) = 1;

        % Q matrix discretization-2 order
        Q_basic = G*qdt*G';
        M1 = Q_basic;
        M2 = Q_basic*F'+F*Q_basic;
        Q = dt*M1 + 1/2*dt*dt*M2;

        % PHIM matrix discretization-2 order
        I = eye(9, 9);
        PHIM = I + dt*F + 1/2*dt*dt*F*F;

        % predict
        x = PHIM*x;
        P = PHIM*P*PHIM' + Q;

        % update from acc
        H = zeros(3, 9);
        H(1, 2) = G_vector(3);
        H(2, 1) = -G_vector(3);
        H(1, 7) = Cbn(1, 1);
        H(1, 8) = Cbn(1, 2);
        H(1, 9) = Cbn(1, 3);
        H(2, 7) = Cbn(2, 1);
        H(2, 8) = Cbn(2, 2);
        H(2, 9) = Cbn(2, 3);
        H(3, 7) = Cbn(3, 1);
        H(3, 8) = Cbn(3, 2);
        H(3, 9) = Cbn(3, 3);

        R = eye(3, 3);
        R(1, 1) = 0.5^2;
        R(2, 2) = 0.5^2;
        R(3, 3) = 0.5^2;
        
        g_estimate = Cbn*(acc_bias - Acc');
        Z = G_vector - g_estimate;
        K = P*H'*((H*P*H'+R)^-1);
        x = x + K*(Z - H*x);
        P = (I - K*H)*P;

        [deltaCbn] = euler2dcm (x(3), x(2), x(1)); % (I+P)Cbn
        Cbn = deltaCbn*Cbn;
        
        % check mag condition
        % calibrated mag data
        Mag_tr = Cali_W_inv*(Mag - Cali_V)';
        mag_update = 1;
        mag_residual(sensor_count, :) = Mag_vector - Cbn*Mag_tr;
        if abs(mag_residual(sensor_count, 1)) > 5
            mag_update = 0;
        end
        
        % update from mag
        if mag_update
        H = zeros(3, 9);
        H(1, 2) = Mag_vector(3);
        H(2, 1) = -Mag_vector(3);
        H(2, 3) = Mag_vector(1);
        H(3, 2) = -Mag_vector(1);

        R = eye(3, 3);
        R(1, 1) = 20^2;
        R(2, 2) = 20^2;
        R(3, 3) = 20^2;

        mag_estimate = Cbn*Mag_tr;
        Z = Mag_vector - mag_estimate;

        K = P*H'*((H*P*H'+R)^-1);
        x = x + K*(Z - H*x);
        P = (I - K*H)*P;
    
        [deltCbn] = euler2dcm (x(3), x(2), x(1)); % (I+P)Cbn
        Cbn = deltCbn*Cbn;
        end
        acc_bias = acc_bias + x(7:9);
        gyro_bias = gyro_bias + x(4:6);
        x(1:9) = 0;

        [yaw, pitch, roll] = dcm2euler(Cbn);
        q = euler2q(yaw, pitch, roll);
        q = q_norm(q);
        heading_fusion(sensor_count) = yaw;
        
        %% step detection
        if step_start_flag ~= 0
            acc_data = Acc;
            acc_det = sqrt(sum(acc_data.^2));
            % smooth acc data
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
        end
        
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
                 if (acc_det_mean - Ge) > 0.8
                     threshold = 0.15 * acc_det_mean;
                 else
                     acc_det_mean = Ge;
                     threshold = 0.08 * acc_det_mean;
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
                    step_detect_flag = 1;
                    step_timing(step_count) = acc_time_tag_array(acc_array_index-1);
                    step_freq = 1 / delta_step_time;
                    delta_step_time = 0;
                else
                    if delta_step_time > 3 % 3s has no step for a long time, need to reset. 
                       delta_step_time = 0;
                    end
                end

                pre_slop = slop;
                pre_step_det = step_det;
        end
        
        %% PDR estimate latitude and longitude
        if step_detect_flag ~= 0
            step_detect_flag = 0;
            step_heading = yaw;
            % calculate the step trajectory
            RM = RE * (1 - esqu) / ((1 - esqu * sin(pdr_latitude) * sin(pdr_latitude))^1.5);
            RN = RE / (sqrt(1 - esqu * sin(pdr_latitude) * sin(pdr_latitude)));

            pdr_latitude = pdr_latitude + step_length*cos(step_heading)/RM;
            pdr_longitude = pdr_longitude + step_length*sin(step_heading)/(RN*cos(pdr_latitude));
            pdr_altitude = gnss_altitude;

            % store all step trajectory
            pdr_latitude_array(step_count) = pdr_latitude*180/pi;
            pdr_longitude_array(step_count) = pdr_longitude*180/pi;
            pdr_altitude_array(step_count) = pdr_altitude;
        end
    end

end

%% convert pdr trajectory to KML file
gps2kml('pdrOut',step_timing,pdr_latitude_array,pdr_longitude_array,pdr_altitude_array,'o-b','MarkerSize',10,'LineWidth',3);

figure;
plot(heading_fusion*180/pi, 'r');
hold on;
plot(heading_pure*180/pi, 'b');
plot(heading_gnss*180/pi, 'g');
title('heading');
legend('fusion', 'pure', 'gnss');

figure;
plot(mag_residual(:, 1), 'r');
hold on;
plot(mag_residual(:, 2), 'g');
plot(mag_residual(:, 3), 'b');
title('mag residual');
legend('x', 'y', 'z');
