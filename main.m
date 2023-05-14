% Clear workspace and command window
clc;
clear;

% Generate the DTMF sequence for phone number 01014227773
% Followed by a silence (guard) interval of 20 ms
number = '01014227773';
fs = 8000; Ts = 1/fs;
x = [];
for i = 1:length(number)
    symbol = number(i);
    dtmf_signal = sym2TT(symbol);
    x = [x, dtmf_signal, zeros(1, fs*0.02)];
end

% Add additive white Gaussian noise with variance 0.1
SNR = 10*log10(var(x)/0.1); % Calculate the required SNR for a variance of 0.1
y = awgn(x, SNR, 'measured'); % Add noise to the signal x

% Write audio file
filename = 'y(t).wav';
myaudiowrite(filename, y, fs,"BitsPerSample", 16);

% Load audio file
[y,fs] = audioread('y(t).wav');

% Plot the signal y(t)
t = (0:1:length(y)-1)/fs;
figure;
plot(t, y);
xlabel('time (s)');
ylabel('y(t)');
title('DTMF Signal with AWGN');

% Compute the signal in the frequency domain
Y = fft(y);
f = linspace(-fs/2,fs/2,length(Y));
Y_Shifted = fftshift(Y);

% Plot the signal in the frequency domain
figure;
plot(f,20*log10(abs(Y_Shifted)))
xlim([600,1700])
title('Spectrum Y(f) of the signal y(t)')
xlabel('f (Hz)')
ylabel('Magnitude (dB)')

% Define the parameters
overlap = 0.5; % 50% overlap
fft_size = 2^14; % FFT size
window_sizes = [16, 64, 256, 1024, 4096]; % Time-domain window sizes

% Create spectrograms using rectangular time-domain window
for i = 1:length(window_sizes)
    window = rectwin(window_sizes(i));
    [S,F,T] = spectrogram(y,window,overlap*length(window),fft_size,fs,'yaxis');
    figure;
    imagesc(T,F,20*log10(abs(S)))
    set(gca,'YDir','normal')
    title(sprintf('Spectrogram (Rectangular Window) with Window Size = %d',window_sizes(i)))
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    colorbar
end

% Worst time-domain resolution: t4096
% Worst frequency-domain resolution: t16
% Optimal trade-off: t256

% Create spectrograms using blackman time-domain window
for i = 1:length(window_sizes)
    window = blackman(window_sizes(i));
    [S,F,T] = spectrogram(y,window,overlap*length(window),fft_size,fs,'yaxis');
    figure;
    imagesc(T,F,20*log10(abs(S)))
    set(gca,'YDir','normal')
    title(sprintf('Spectrogram (Blackman Window) with Window Size = %d',window_sizes(i)))
    xlabel('Time (s)')
    ylabel('Frequency (Hz)')
    colorbar
end

% Blackman window provides better frequency resolution than rectangular window

% Decode DTMF symbols using Goertzel algorithm
symbol_length = 0.1; % Duration of each DTMF symbol
guard_length = 0.02; % Guard interval
start_idx = 1;
end_idx = 0;
symbols_detected = [];
while start_idx <= length(y) - fs*symbol_length
    % Extract the symbol signal
    end_idx = start_idx + fs*symbol_length - 1;
    symbol_signal = y(start_idx:end_idx); % The segment extracted that corresponds to a single DTMF symbol

    % Decode the symbol using the Goertzel algorithm
    DTMF_frequencies = [697, 770, 852, 941, 1209, 1336, 1477, 1633];
    DTMF_map = [1, 2, 3; 4, 5, 6; 7, 8, 9; -1, 0, -2];
    threshold = 10e3;
    detected_frequencies = [];
    for i = 1:length(DTMF_frequencies)
        freq = DTMF_frequencies(i);
        k = freq*symbol_length/fs; % Compute the index of the desired frequency bin
        w = 2*pi*k/symbol_length;
        coeff = 2*cos(w);
        Q2 = 0;
        Q1 = 0;
        for j = 1:length(symbol_signal)
            Q0 = symbol_signal(j) + coeff*Q1 - Q2;
            Q2 = Q1;
            Q1 = Q0;
        end
        power = abs(Q1)^2 + abs(Q2)^2 - 0.5*coeff*Q1*Q2;
        % Power calculated for the frequency bin with index (k) of the desired frequency (freq)
        if power > threshold
            detected_frequencies = [detected_frequencies, freq];
        end
    end
    if length(detected_frequencies) == 2 % Only 2 frequencies detected, referring to a specific DTMF symbol
        row = mod(find(DTMF_frequencies == detected_frequencies(1)),4);
        col = mod(find(DTMF_frequencies == detected_frequencies(2)),4);
        if row == 0
            row = 4;
        end
        if col == 0
            col = 4;
        end
        symbol = DTMF_map(row, col);
        symbols_detected = [symbols_detected, symbol];
    end

    start_idx = end_idx + fs*guard_length; % start_idx is a 'symbol_length + guard_length' over the current one
end

disp('Symbols detected using Goertzel algorithm:');
disp(symbols_detected);