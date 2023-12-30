%Read the audio file
[y, Fs] = audioread('eric.wav'); 
t1 = linspace(0, length(y)/Fs , length(y)) ;
figure
subplot(2,1,1)
plot (t1 ,y ,'b');
xlabel('time(s)') ;
title('Waveform of original signal');
% Find the spectrum of the signal
N = length(y); % Number of samples
Y = fft(y); % Fast Fourier transform
f = (-N/2:N/2-1)*Fs/N; % Frequency vector
Y = fftshift(Y); % Shift zero frequency to the center
% Plot the spectrum
subplot(2,1,2)
plot(f, abs(Y)/Fs,'b')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of the original signal')

% Remove all frequencies greater than 4 kHz
fc = 4000; % Cutoff frequency
Y(abs(f)  > fc) = 0; % Set high frequencies to zero
% Obtain the filtered signal in time and frequency domains
yf = ifft(ifftshift(Y)); % Inverse fast Fourier transform
Yf = fft(yf); % Fast Fourier transform of the filtered signal
Yf = fftshift(Yf); % Shift zero frequency to the center
% Plot the spectrum of the filtered signal
Yf_T = real(ifft(ifftshift(Yf))) ;
figure
subplot(2,1,1)
plot(t1, Yf_T,'b');
title('Waveform of filtered signal') ;
subplot(2,1,2)
plot(f, abs(Yf)/Fs,'b')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of the filtered signal')
% Play the filtered signal 
%Yf_T = Yf_T / max(abs(Yf_T)); %normalize the amplitude of your modulating signal,
player = audioplayer(Yf_T, Fs);
%play(player);
%pause(player.TotalSamples / Fs);

% NBFM Modulation
Fc = 100000; % carrier frequency
Fs2 = 5 * Fc; % sampling frequency
Kf = 0.1; % frequency Sensitivity
A = 1; % amplitude

% resampling the msg
msg_resampled = resample(Yf_T , Fs2 , Fs);

% generating linspaces
t2 = linspace (0 ,  (length(msg_resampled)) / Fs2 , length(msg_resampled)) ;
f2 = linspace( (-Fs2)/2 , Fs2/2 , length(msg_resampled) ) ;


NBFM_signal_t = A * cos(2*pi*Fc* t2)' - A * 2*pi*Kf * sin(2*pi*Fc* t2)' .* cumsum(msg_resampled);
NBFM_signal_freq = fftshift ( fft (NBFM_signal_t)) ;

% plots
figure
subplot(2,1,1);
plot ( t2,NBFM_signal_t);
title('Waveform of FM modulation msg');
subplot(2,1,2);
plot(f2, abs(NBFM_signal_freq)/Fs2);
title('Spectrum');

%demodulation

env = abs(hilbert(NBFM_signal_t));

% differentiator
dem_msg = [0; diff(env)];
% DC blocker
dem_msg = dem_msg - mean(dem_msg);
%resample
dem_msg = 1.5*resample(dem_msg, Fs, Fs2);
dem_msg = dem_msg(1:length(y));

%plot
figure
subplot(2,1,1)
plot(t1, dem_msg); 
xlabel('time(s)') ;
title('Demodulated Message') ;
ylim([-0.4 0.2]);
Y2 = fft(dem_msg);
Y2 = fftshift(Y2); % Shift zero frequency to the center
% Plot the spectrum
subplot(2,1,2)
plot(f, abs(Y2)/Fs,'b')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
title('Spectrum of the demodualted signal')
%playing final sound
sound(dem_msg, Fs);


