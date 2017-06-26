function[frequency, amplitude] = fft_process(data, sample_rate, sample_time)

addpath('.\functions\');

%% smooth the data, get average value every 5 points, and the sample_rate decrease 5 times
count = 1;
ave_num = 5;
sum = 0;

for i = 1:length(data)
   if mod(i, ave_num) == 0
       filtered_data(count) = (sum + data(i)) / ave_num;
       count = count + 1;
       sum = 0;
   else
       sum = sum + data(i);
   end
end

%% Remove constant signal
det_mean = mean(filtered_data);
filtered_data = filtered_data - det_mean;

%% fft process
Fs = sample_rate / ave_num;   % Sampling frequency
T = 1/Fs;
L = length(filtered_data);    % Length of signal
t = (0:L-1)*T;                % Time vector
y = filtered_data;
NFFT = 2^nextpow2(L);         % Next power of 2 from length of y
Y = fft(y,NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);
[value, index] = max(2*abs(Y));

amplitude = value;
frequency = f(index);


%% plot figure for debug
if 0
    figure;
    plot(t(1:Fs*sample_time),y(1:Fs*sample_time))
    title('Filtered Signal in 3 Second')
    xlabel('time (seconds)')

    % Plot single-sided amplitude spectrum.
    figure;
    plot(f,2*abs(Y(1:NFFT/2+1))) 
    title('Single-Sided Amplitude Spectrum of Filetered Signal')
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
end