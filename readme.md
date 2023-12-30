# FM modulation

This MATLAB script demonstrates the process of frequency modulation (FM) and its demodulation. Let's break down the code:

1. **Read the Audio File:**
   - Reads an audio file named 'eric.wav'.
   - Plots the waveform and spectrum of the original signal.

```matlab
[y, Fs] = audioread('eric.wav');
t1 = linspace(0, length(y)/Fs, length(y));
figure
subplot(2,1,1)
plot(t1, y, 'b');
xlabel('time(s)');
title('Waveform of original signal');

N = length(y); % Number of samples
Y = fft(y); % Fast Fourier transform
f = (-N/2:N/2-1)*Fs/N; % Frequency vector
Y = fftshift(Y); % Shift zero frequency to the center
subplot(2,1,2)
plot(f, abs(Y)/Fs, 'b')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of the original signal')
```

2. **Filtering:**
   - Removes frequencies greater than 4 kHz.
   - Plots the waveform and spectrum of the filtered signal.

```matlab
fc = 4000; % Cutoff frequency
Y(abs(f) > fc) = 0; % Set high frequencies to zero
yf = ifft(ifftshift(Y)); % Inverse fast Fourier transform
Yf = fft(yf); % Fast Fourier transform of the filtered signal
Yf = fftshift(Yf); % Shift zero frequency to the center
Yf_T = real(ifft(ifftshift(Yf)));
figure
subplot(2,1,1)
plot(t1, Yf_T, 'b');
title('Waveform of filtered signal');
subplot(2,1,2)
plot(f, abs(Yf)/Fs, 'b')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of the filtered signal')
```

3. **NBFM Modulation:**
   - Modulates the filtered signal using Narrowband Frequency Modulation (NBFM).
   - Plots the waveform and spectrum of the modulated signal.

```matlab
Fc = 100000; % Carrier frequency
Fs2 = 5 * Fc; % Sampling frequency
Kf = 0.1; % Frequency Sensitivity
A = 1; % Amplitude

msg_resampled = resample(Yf_T, Fs2, Fs);
t2 = linspace(0, length(msg_resampled)/Fs2, length(msg_resampled));
f2 = linspace((-Fs2)/2, Fs2/2, length(msg_resampled));

NBFM_signal_t = A * cos(2*pi*Fc*t2)' - A * 2*pi*Kf * sin(2*pi*Fc*t2)' .* cumsum(msg_resampled);
NBFM_signal_freq = fftshift(fft(NBFM_signal_t));

figure
subplot(2,1,1);
plot(t2, NBFM_signal_t);
title('Waveform of FM modulation msg');
subplot(2,1,2);
plot(f2, abs(NBFM_signal_freq)/Fs2);
title('Spectrum');
```

4. **Demodulation:**
   - Demodulates the NBFM signal back to the original signal.
   - Plots the waveform and spectrum of the demodulated signal.

```matlab
env = abs(hilbert(NBFM_signal_t));
dem_msg = [0; diff(env)];
dem_msg = dem_msg - mean(dem_msg);
dem_msg = 1.5 * resample(dem_msg, Fs, Fs2);
dem_msg = dem_msg(1:length(y));

figure
subplot(2,1,1)
plot(t1, dem_msg); 
xlabel('time(s)');
title('Demodulated Message');
ylim([-0.4 0.2]);
Y2 = fft(dem_msg);
Y2 = fftshift(Y2);
subplot(2,1,2)
plot(f, abs(Y2)/Fs,'b')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of the demodulated signal')
```

5. **Playing the Final Sound:**
   - Plays the demodulated message.

```matlab
sound(dem_msg, Fs);
```

Note: The modulation and demodulation processes involve resampling and filtering. The script then visualizes the results in both time and frequency domains. Finally, it plays the demodulated sound.