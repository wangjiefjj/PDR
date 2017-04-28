clear all;
clc;

%% load acc data
data = load('.\mrvlAcc_walking_50hz.txt');
acc_data = data(:, 3:5);
acc_det = sqrt(acc_data(:,1).^2 + acc_data(:,2).^2 + acc_data(:,3).^2);
sample_rate = 50;
ave_num = 5;

%% smooth the data, get average value every 5 points, and the sample_rate decrease 5 times
count = 1;
sum = 0;
len = floor(length(acc_data) / ave_num);
acc_det_filtered = zeros(len, 1);

for i = 1:length(acc_data)
   if mod(i, ave_num) == 0
       acc_det_filtered(count) = (sum + acc_det(i)) / ave_num;
       count = count + 1;
       sum = 0;
   else
       sum = sum + acc_det(i);
   end
end

%% plot the filter acc det data for debug
Fs = sample_rate / ave_num;   % Sampling frequency
T = 1/Fs;
L = length(acc_det_filtered); % Length of signal
t = (0:L-1)*T;                % Time vector
y = acc_det_filtered;

figure;
start_time = 100*Fs + 1;
sample_time = 10;             % Display "sample_time" length data
plot(t(start_time:start_time + Fs*sample_time),y(start_time:start_time + Fs*sample_time))
title('Filtered Acc Det in 3 Second')
xlabel('Time (seconds)')
ylabel('Y(t)')

%% start step detection (looking for wave summit)
window_length = 50;
acc_det_mean = 0;
Gravity = 9.8;
pre_step_det = Gravity;
pre_slop = 0;
step_lock_count = round(0.3*Fs);
delta_step_count = 0;
step_count = 0;

% signal_end = length(acc_det_filtered);
signal_start = start_time;
signal_end = start_time + Fs*sample_time;
for i = signal_start:signal_end
    step_det = acc_det_filtered(i);
    % average calculate and the window length is 50
    if i < window_length
        acc_det_mean = (acc_det_mean*(i-1) + step_det)/i;
    else
        acc_det_mean = (acc_det_mean*(window_length-1) + step_det)/window_length;
    end

    % determin the threshold
    if abs(acc_det_mean - Gravity) > 0.8
        threshold = 0.3 * acc_det_mean;
    else
        acc_det_mean = Gravity;
        threshold = 0.1 * acc_det_mean;
    end
    
    % slop calculate
    if step_det - pre_step_det > 0
        slop = 1;
    else
        slop = -1;
    end
    
    % step detect (wave summit & delta time & step det)
    if ( pre_slop > 0 && slop < 0  && delta_step_count >= step_lock_count && pre_step_det - acc_det_mean > threshold)
        step_count = step_count + 1;
        step_timing(step_count) = t(i-1);
        delta_step_count = 0;
    else
        delta_step_count = delta_step_count + 1;
        if delta_step_count > Fs*10 % 10s
           % has no step for a long time, need to reset the step detection module. 
           delta_step_count = 0;
        end
    end
    
    pre_slop = slop;
    pre_step_det = step_det; 
    
end

fprintf('step count = %d', step_count);
















