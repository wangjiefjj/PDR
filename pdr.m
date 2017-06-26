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
        Mag(count, :) = data(i, 9:11);
        Acc(count, :) = data(i, 3:5);
        Gyro(count, :) = data(i, 6:8);
        count = count + 1;
    end
end

% disp gyro data
if 0
    figure;
    plot(Gyro(:, 1), 'r');
%     hold on;
%     plot(Gyro(:, 2), 'g');
%     plot(Gyro(:, 3), 'b');
    legend('X', 'Y', 'Z');
    title('gyro data');
end

% mag_calibrated = inv(W)*(mag_raw - V)
N = SampleRate*MagCalibrationTime;
[Cali_B, Cali_V, Cali_W_inv, Cali_Error] = cali7eig(Mag(1:N, :));
if Cali_Error < 0.15
    mag_calibration = 1;
else
    mag_calibration = 0;
    disp('mag calibration can not be execuated');
end

%% configure

SampleRate = 50;    % Hz
AlignmentTime = 2;  % second

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
mag_alignment = Mag(M:M+N, :);

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

%% sensor fusion variable

Ge = 9.80665;
G_vector = [0, 0, Ge]';
yaw = zeros(N, 1);
pitch = zeros(N, 1);
roll = zeros(N, 1);
gyro_bias = zeros(3, 1);
acc_bias = zeros(3, 1);
mag_jamming = zeros(3, 1);
Cnb = eye(3, 3);
Cbn = eye(3, 3);

StateNum = 12;
MeasNumG = 3;
MeasNumM = 3;
x = zeros(StateNum, 1); % roll, pitch, yaw, gyro_bias_x, gyro_bias_y, gyro_bias_z, acc_bias_x, acc_bias_y, acc_bias_z, mag_jamming_x, mag_jamming_y, mag_jamming_z
Corr_time_gyro = 1;
Corr_time_acc = 1;
Corr_time_mag = 10;
sigma_Win = 1.0e-6;
sigma_gyro_1 = (2*pi/180/3600)^2;   % Markov process
sigma_gyro_2 = (10*pi/180/3600)^2;   % Random walk
sigma_acc = ((5.0e-2)*Ge)^2;  % Markov process
sigma_mag = 0.1*0.1;  % Markov process

error_phim_e = 1.0*pi/180;
error_phim_n = 1.0*pi/180;
error_phim_u = 1.0*pi/180;
error_gyro = 1000*pi/180/3600;
error_acc = 0.3;
error_mag_jamming = 5;
P = zeros(StateNum, StateNum);
P(1, 1) = error_phim_e^2;
P(2, 2) = error_phim_n^2;
P(3, 3) = error_phim_u^2;
P(4, 4) = error_gyro^2;
P(5, 5) = error_gyro^2;
P(6, 6) = error_gyro^2;
P(7, 7) = error_acc^2;
P(8, 8) = error_acc^2;
P(9, 9) = error_acc^2;
P(10, 10) = error_mag_jamming^2;
P(11, 11) = error_mag_jamming^2;
P(12, 12) = error_mag_jamming^2;

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
step_length = 0.7;

%% gnss variable
RE = 6378137.0;
esqu = 0.00669437999013;

%% pdr variable
step_start_flag = 0;

%% step length variable
step_adjust_window = 10;
step_adjust_count = 0;

%% main loop start
gnss_count = 1;
sensor_count = 0;
gyro_smooth_count = 1;
gyro_smooth_number = 10;
gyro_smooth_array = zeros(gyro_smooth_number, 3);
heading_smooth_count = 1;
heading_smooth_number = 5;
heading_smooth_array = zeros(heading_smooth_number, 1);

for i = M:length(data)
    type = data(i, 1);
    time_tag = data(i, 2);
    % GNSS
    if type == GNSS
        gnss_latitude = data(i, 3);
        gnss_longitude = data(i, 4);
        gnss_altitude = data(i, 5);
        gnss_heading(gnss_count) = data(i, 9);
        if gnss_count == 1
            q = euler2q(yaw_initial, pitch_initial, roll_initial);
            qlpf = q;
            qlpf_6D = q;
            sensor_heading_9D(1) = yaw_initial;
            sensor_heading_6D(1) = yaw_initial;
            mag_inclination(1) = geo_inclination;
            mag_inclination_lpf = geo_inclination;
            % initial the latitude and longitude
            step_start_flag = 1;
            pdr_latitude = gnss_latitude;
            pdr_longitude = gnss_longitude;
            pdr_altitude = gnss_altitude;
        end
        if sensor_count > 0
            % sensor fusion
            dt = 1/SampleRate*sensor_count;
            Acc = mean(acc_array);
            Mag = mean(mag_array);
            
           %% mag + acc heading estimation
            Cnb_6D = ecompass_ned(-Acc, Mag);
            Cbn_6D = Cnb_6D';
            [yaw_6D, pitch_6D, roll_6D] = dcm2euler(Cbn_6D);
            q_6D = euler2q(yaw_6D, pitch_6D, roll_6D);
            q_6D = q_norm(q_6D);
            % low pass filter for quaternion
            % set low pass filter constant with maximum value 1.0 (all pass) decreasing to 0.0 (increasing low pass)
            lpf_time = 5;   % time constant (second)
            flpf = dt/lpf_time;
            deltaq = qconjgAxB(qlpf_6D, q_6D);
            if deltaq(1) < 0
                deltaq = -deltaq;
            end
            ftemp = flpf + (1-flpf)*(1-deltaq(1));
            deltaq = ftemp*deltaq;
            deltaq(1) = sqrt(1 - norm(deltaq(2:4))^2);
            qlpf_6D = qAxB(qlpf_6D, deltaq);
            qlpf_6D = q_norm(qlpf_6D);
            % convert to euler
            Cbn_6D = q2dcm(qlpf_6D);
            [sensor_heading_6D(gnss_count), ~, ~] = dcm2euler(Cbn_6D);
            
           %% gyro + acc + mag heading estimation 
            if 1
            % kalman sensor fusion
            Cbn = q2dcm(q);

            Fg = diag([-1/Corr_time_gyro, -1/Corr_time_gyro, -1/Corr_time_gyro]);
            Fa = diag([-1/Corr_time_acc, -1/Corr_time_acc, -1/Corr_time_acc]);
            Fm = diag([-1/Corr_time_mag, -1/Corr_time_mag, -1/Corr_time_mag]);
            Fnull = zeros(3, 3);

            F = [Fnull, -Cbn, Fnull, Fnull;
                 Fnull, Fg,   Fnull, Fnull;
                 Fnull, Fnull,Fa,    Fnull;
                 Fnull, Fnull,Fnull, Fm  ];

            qdt = diag([sigma_Win, sigma_Win, sigma_Win, sigma_gyro_1, sigma_gyro_1, sigma_gyro_1, sigma_gyro_2,...
                        sigma_gyro_2, sigma_gyro_2, sigma_acc, sigma_acc, sigma_acc, sigma_mag, sigma_mag, sigma_mag]);
            I = eye(3, 3);
            Gnull = zeros(3, 3);
            G = [-Cbn,    Gnull,  -Cbn,  Gnull,  Gnull;
                  Gnull,      I,  Gnull, Gnull,  Gnull;
                  Gnull,  Gnull,  Gnull,     I,  Gnull;
                  Gnull,  Gnull,  Gnull, Gnull,      I];
            % Q matrix discretization-2 order
            Q_basic = G*qdt*G';
            M1 = Q_basic;
            M2 = Q_basic*F'+F*Q_basic;
            Q = dt*M1 + 1/2*dt*dt*M2;

            % PHIM matrix discretization-2 order
            I = eye(StateNum, StateNum);
            PHIM = I + dt*F + 1/2*dt*dt*F*F;

            % predict
            x = PHIM*x;
            P = PHIM*P*PHIM' + Q;
            
            % update from acc
            H = zeros(3, StateNum);
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

            R = eye(MeasNumG, MeasNumG);
            R(1, 1) = 5^2;
            R(2, 2) = 5^2;
            R(3, 3) = 5^2;

            acc_liner_b = zeros(3, 1);
            g_estimate = Cbn*(acc_bias - Acc' + acc_liner_b);
            Z = G_vector - g_estimate;
            K = P*H'*((H*P*H'+R)^-1);
            x = x + K*(Z - H*x);
            P = (I - K*H)*P;

            [deltaCbn] = euler2dcm (x(3), x(2), x(1)); % (I+P)Cbn
            Cbn = deltaCbn*Cbn;
            
            % update from mag
            H = zeros(3, StateNum);
            H(1, 2) = Mag_vector(3);
            H(2, 1) = -Mag_vector(3);
            H(2, 3) = Mag_vector(1);
            H(3, 2) = -Mag_vector(1);
            H(1, 10) = -Cbn(1, 1);
            H(1, 11) = -Cbn(1, 2);
            H(1, 12) = -Cbn(1, 3);
            H(2, 10) = -Cbn(2, 1);
            H(2, 11) = -Cbn(2, 2);
            H(2, 12) = -Cbn(2, 3);
            H(3, 10) = -Cbn(3, 1);
            H(3, 11) = -Cbn(3, 2);
            H(3, 12) = -Cbn(3, 3);

            R = eye(3, 3);
            R(1, 1) = 1^2;
            R(2, 2) = 1^2;
            R(3, 3) = 1^2;

            mag_estimate = Cbn*Mag';
            Z = Mag_vector - mag_estimate;

            K = P*H'*((H*P*H'+R)^-1);
            x = x + K*(Z - H*x);
            
            % mag jamming check
            mag_jamming_b = x(10:12);
            if norm(mag_jamming_b) > 2*Cali_B^2
                disp('mag jamming occur');
            else
                P = (I - K*H)*P;
                [deltaCbn] = euler2dcm (x(3), x(2), x(1)); % (I+P)Cbn
                Cbn = deltaCbn*Cbn;
                [yaw(i), pitch(i), roll(i)] = dcm2euler(Cbn);   % gyro + acc + mag estimation
                q = euler2q(yaw(i), pitch(i), roll(i));
                q = q_norm(q);

                % mag vector correction
                mag_jamming_n = Cbn*mag_jamming_b;
                mag_vector_temp = Mag_vector - mag_jamming_n;
                geoB_temp = sqrt(mag_vector_temp(1)^2 + mag_vector_temp(3)^2);
                sindelta = mag_vector_temp(3)/geoB_temp;
                cosdelta = mag_vector_temp(1)/geoB_temp;
                SINDELTAMAX = 0.9063078; % sin of max +ve geomagnetic inclination angle: here 65.0 deg
                COSDELTAMAX = 0.4226183; % cos of max +ve geomagnetic inclination angle: here 65.0 deg
                if sindelta > SINDELTAMAX || sindelta < -SINDELTAMAX
                    sindelta = SINDELTAMAX;
                    cosdelta = COSDELTAMAX;
                end
                geo_inclination = asin(sindelta);
                % low pass filter
                if 1
                    lpf_time = 5;   % time constant (second)
                    flpf = dt/lpf_time;
                    mag_inclination_lpf = mag_inclination_lpf + flpf * (geo_inclination - mag_inclination_lpf);
                    geo_inclination = mag_inclination_lpf;
                end
                mag_inclination(gnss_count) = geo_inclination;
                if mag_calibration == 0
                    Cali_B = geoB_temp;
                end
                Mag_vector = [Cali_B*cos(geo_inclination), 0, Cali_B*sin(geo_inclination)]';
            end

            [deltaCbn] = euler2dcm (x(3), x(2), x(1)); % (I+P)Cbn
            Cbn = deltaCbn*Cbn;
            [yaw, pitch, roll] = dcm2euler(Cbn);   % gyro + acc + mag estimation
            q = euler2q(yaw, pitch, roll);
            q = q_norm(q);
            
            % feedback
            acc_bias = acc_bias + x(7:9);
            gyro_bias = gyro_bias + x(4:6);
            x(1:StateNum) = 0;
            end
            
            if 0
            % low pass filter for quaternion
            % set low pass filter constant with maximum value 1.0 (all pass) decreasing to 0.0 (increasing low pass)
            lpf_time = 1;   % time constant (second)
            flpf = dt/lpf_time;
            deltaq = qconjgAxB(qlpf, q);
            if deltaq(1) < 0
                deltaq = -deltaq;
            end
            ftemp = flpf + (1-flpf)*(1-deltaq(1));
            deltaq = ftemp*deltaq;
            deltaq(1) = sqrt(1 - norm(deltaq(2:4))^2);
            qlpf = qAxB(qlpf, deltaq);
            qlpf = q_norm(qlpf);
            % convert to euler
            Cbn = q2dcm(qlpf);
            [yaw, ~, ~] = dcm2euler(Cbn);
            q = qlpf;
            end
            
            % heading smooth
            if heading_smooth_count > heading_smooth_number
                for j = 1 : heading_smooth_number - 1
                    heading_smooth_array(j) = heading_smooth_array(j+1);
                end
                heading_smooth_array(heading_smooth_number) = yaw;
            else
                heading_smooth_array(heading_smooth_count) = yaw;
            end
            slot_number = min(heading_smooth_count, heading_smooth_number);
            temp_array = zeros(slot_number, 1);
            % check discontinuous
            temp_array(1) = heading_smooth_array(1);
            for j = 2:slot_number
                temp_array(j) = heading_smooth_array(j);
                if temp_array(j) - temp_array(j-1) > 180/180*pi
                    temp_array(j) = temp_array(j) - 2*pi;
                end
                if temp_array(j) - temp_array(j-1) < -180/180*pi
                    temp_array(j) = temp_array(j) + 2*pi;
                end
            end
            yaw = mean(temp_array(1:slot_number));
            if yaw > pi
                yaw = yaw - 2*pi;
            end
            if yaw < -pi
                yaw = yaw + 2*pi;
            end
            heading_smooth_count = heading_smooth_count + 1;
            
            % restore the heading result
            sensor_heading_9D(gnss_count) = yaw;
        else
            if gnss_count ~= 1
                % lost sensor data
                sensor_heading_9D(gnss_count) = sensor_heading_9D(gnss_count - 1);
                sensor_heading_6D(gnss_count) = sensor_heading_6D(gnss_count - 1);
                mag_inclination(gnss_count) = mag_inclination(gnss_count - 1);
            end
        end
        
        gnss_count = gnss_count + 1;
        sensor_count = 0;
    % SENSOR
    elseif type == SENSOR
        sensor_count = sensor_count + 1;
        % store sensor data into pool
        Mag = data(i, 9:11);
        Acc = data(i, 3:5);
        Gyro = data(i, 6:8);
        % calibrated mag data
        Mag = Cali_W_inv*(Mag - Cali_V)';
        mag_array(sensor_count, :) = Mag;
        acc_array(sensor_count, :) = Acc;
        % smooth gyro data
        if gyro_smooth_count > gyro_smooth_number
            for j = 1 : gyro_smooth_number - 1
                gyro_smooth_array(j, :) = gyro_smooth_array(j+1, :);
            end
            gyro_smooth_array(gyro_smooth_number, :) = Gyro; 
        else
            gyro_smooth_array(gyro_smooth_count, :) = Gyro;
        end
        slot_number = min(gyro_smooth_count, gyro_smooth_number);
        Gyro(1) = mean(gyro_smooth_array(1:slot_number, 1));
        Gyro(2) = mean(gyro_smooth_array(1:slot_number, 2));
        Gyro(3) = mean(gyro_smooth_array(1:slot_number, 3));
        gyro_smooth_count = gyro_smooth_count + 1;
        % strapdown mechanization
        Cbn = q2dcm(q);
        Cnb = Cbn';

        Wepp = zeros(3, 1); % no latitude information in computing latitude and longitude rate
        Wiep = zeros(3, 1); % no latitude information in computing earth rate in the navigation frame
        Wipp = Wiep + Wepp;
        Wipb = Cnb * Wipp;
        Wpbb = Gyro' - gyro_bias - Wipb;
        
        for j = 1:3
            if abs(Wpbb(j)) < 30 / 180 * pi
                Wpbb(j) = 0;
            end
        end

        dq = zeros(4, 1);
        dq(1) = -(Wpbb(1)*q(2) + Wpbb(2)*q(3) + Wpbb(3)*q(4))/2;
        dq(2) = (Wpbb(1)*q(1) + Wpbb(3)*q(3) - Wpbb(2)*q(4))/2;
        dq(3) = (Wpbb(2)*q(1) - Wpbb(3)*q(2) + Wpbb(1)*q(4))/2;
        dq(4) = (Wpbb(3)*q(1) + Wpbb(2)*q(2) - Wpbb(1)*q(3))/2;
        dt = 1/SampleRate;
        q = q + dq*dt;
        q = q_norm(q);
        Cbn = q2dcm(q);
        [yaw, pitch, roll] = dcm2euler(Cbn); % pure gyro estimation
        
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
                 if abs(acc_det_mean - Ge) > 0.8
                     threshold = 0.15 * acc_det_mean;
                 else
                     acc_det_mean = Ge;
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
end

%% plot the filter acc det data for debug
Fs = SampleRate/ave_num;            % About 10Hz
L = length(acc_det_filtered_array); % Length of signal
t = acc_time_tag_array;             % Time vector
y = acc_det_filtered_array;

if 1
figure;
start_time = 0*Fs + 1;
sample_time = 10;                  % Display "sample_time" length data
plot(t, y);
% plot(t(start_time:start_time + Fs*sample_time),y(start_time:start_time + Fs*sample_time))
title('Filtered Acc Det in 10 Second')
xlabel('Time (seconds)')
ylabel('Y(t)')
end

figure;
plot(gnss_heading*180/pi, 'r');
hold on;
plot(sensor_heading_9D*180/pi, 'g');
plot(sensor_heading_6D*180/pi, 'b');
legend('gnss heading', '9D heaing', '6D heading');
ylabel('angle(degree)');
xlabel('time(second)');
title(['heading comparison, 9D heading std = ', num2str(std(sensor_heading_9D*180/pi)),'deg, ', '6D heading std = ', num2str(std(sensor_heading_6D*180/pi)),'deg']);

figure;
plot(mag_inclination*180/pi);
ylabel('angle(degree)');
xlabel('time(second)');
title(['geo inclination, std = ', num2str(std(mag_inclination*180/pi)),'deg']);

%% convert pdr trajectory to KML file
gps2kml('pdrOut',step_timing,pdr_latitude_array,pdr_longitude_array,pdr_altitude_array,'o-b','MarkerSize',10,'LineWidth',3);











