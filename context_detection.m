addpath('.\functions\');

clear all;
clc;

STATIC = 0;
WALKING = 1;
DRIVING = 2;
UNKNOWN = -1;
previous_mode = UNKNOWN;

%% load acc data and gnss data
data = load('.\data\context_data\mrvlAcc_roadtest2_50hz.txt');
acc_data = data(:, 3:5);

%% process sample data every 3 second
sample_rate = 50;
sample_time = 3;
signal_length = sample_rate*sample_time;
count = 0;
acc_det_pool_3s = zeros(sample_rate, 1);

result_count = 1;
for i = 1:length(acc_data)
    %% detect static/walk/driving every 3 second
    acc_det = sqrt(acc_data(i,1).^2 + acc_data(i,2).^2 + acc_data(i,3).^2);
    count = count + 1;
    acc_det_pool_3s(count) = acc_det;
    if count == signal_length
        count = 0;
        %% detect static by acc det std
        acc_det_mean = mean(acc_det_pool_3s);
        acc_det_std(result_count) = std(acc_det_pool_3s);
        acc_det_no_bias = acc_det_pool_3s - acc_det_mean;
        acc_det_energy(result_count) = sum(acc_det_no_bias.^2)/signal_length;
        if (acc_det_std(result_count) < 0.25 || acc_det_energy(result_count) < 0.08)
            motion_mode(result_count) = STATIC;
        else
            motion_mode(result_count) = UNKNOWN;
        end
        %% detect walking/driving in frequency domain
        [Frequency, Amplitude] = fft_process(acc_det_pool_3s, sample_rate, sample_time);
        Ft(result_count) = Frequency;
        Amp(result_count) = Amplitude;
        if motion_mode(result_count) ~= STATIC
            if Ft(result_count) > 1.5 && Ft(result_count) < 1.9 && Amp(result_count) > 1.0
                motion_mode(result_count) = WALKING;
            else if Amp(result_count) > 0.2 && Amp(result_count) < 0.6
                    motion_mode(result_count) = DRIVING;
                end
            end
        end
        
        %% include previous information
        if previous_mode == STATIC
            if (acc_det_std(result_count) < 0.9 || acc_det_energy(result_count) < 0.8)
                motion_mode(result_count) = STATIC;
            end
        end

        if previous_mode == WALKING
            if Ft(result_count) > 0.6 && Ft(result_count) < 2.3 || Amp(result_count) > 0.35
                    motion_mode(result_count) = WALKING;
            end
        end

        if previous_mode == DRIVING
            if Amp(result_count) > 0.04 && Amp(result_count) < 2.0
                motion_mode(result_count) = DRIVING;
            end
        end
        previous_mode = motion_mode(result_count);
        result_count = result_count + 1;
    end
     
end

%% Display result
if 0
figure;
plot(Ft, Amp, '.');
title('Single-Sided Amplitude Spectrum of Filetered Signal')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
end

if 1
figure;
plot(acc_det_std, 'r');
hold on;
plot(acc_det_energy, 'b');
legend('acc det std', 'acc det energy');
title('accelerometer data energy and standard deviation');
end

if 1
figure;
plot(motion_mode);
title('motion detect result') 
end

if 1
figure;
plot(Amp);
title('amp')
end

if 0
figure;
plot(Ft)
title('frequency');
end


