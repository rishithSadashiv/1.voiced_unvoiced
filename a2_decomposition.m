% Iterative method for decomposition
% Original decomposition (for plots)


clc;
clear all;
close all;

[d,fs1]=audioread('./arctic_b0399.wav');
disp(length(d));
fs=fs1;
% fs=8000;
% d=resample(d,fs,fs1);
% disp(length(d));

% [d,fs]=audioread('new_wavs/arctic_a0023.wav');
d=d-mean(d);
d=d./max(abs(d));


voiced_index = 5520;%11648;
unvoiced_index = 22336;
ms30 = 0.032*fs;

% Voiced segment (32 ms)
voiced = d(voiced_index:voiced_index+ms30-1);
% LP residual of voiced segment

framesize = 32;
frameshift = 4;
lporder = 10;
res_d=LPres(d,fs,framesize,frameshift,lporder,1);
% res_voiced=LPres(voiced,fs,framesize,frameshift,lporder,1);  % because they have used 4ms shift in paper. Kept frame size same as before.
% res_voiced = res_voiced.*hamming(length(res_voiced));
% [cep_voiced,ym] = rceps(res_voiced);

%%

% figure(1);
% subplot(411);
% plot(voiced);
% subplot(412);
% plot(res_voiced);
% subplot(413);
% plot(res_voiced_w);
% subplot(414);
% plot(cep_voiced_w);

res_voiced = res_d(voiced_index:voiced_index+ms30-1);

N = 1024;
res_X = fft(res_voiced, N);
mag_X = abs(res_X);
phase_X = angle(res_X);
logmag_X = log(mag_X);
cep_x = ifft(logmag_X, N);


logmag_X_dash = fft(cep_x, N);
mag_X_dash = exp(logmag_X_dash);
res_X_dash = mag_X_dash.*exp(1i*phase_X);
res_voiced_dash = real(ifft(res_X_dash, N));

ltf_sample = 0.002*fs; % Low time lifter cutoff 0.002 s
cep_x1 = cep_x;
cep_x1(1:ltf_sample) = 0;
cep_x1(end-ltf_sample+1:end) = 0;
% [pks, locs] = findpeaks(cep_x1);
[vals, locs] = maxk(cep_x1, 2);
% Harmonic part of cepstrun
cep_x1_h = zeros(length(cep_x1),1);
cep_x1_h(locs(1)-1:locs(1)+1) = cep_x1(locs(1)-1:locs(1)+1);
cep_x1_h(2*locs(1)-1:2*locs(1)+1) = cep_x1(2*locs(1)-1:2*locs(1)+1);
% cep_x1_h(3*locs(1)-1:3*locs(1)+1) = cep_x1(3*locs(1)-1:3*locs(1)+1);
cep_x1_h(locs(2)-1:locs(2)+1) = cep_x1(locs(2)-1:locs(2)+1);
cep_x1_h((2*locs(2)-N-2)-1:(2*locs(2)-N-2)+1) = cep_x1((2*locs(2)-N-2)-1:(2*locs(2)-N-2)+1);
% cep_x1_h((3*locs(2)-N-2)-1:(3*locs(2)-N-2)+1) = cep_x1((3*locs(2)-N-2)-1:(3*locs(2)-N-2)+1);
% Noise part of cepstrum
cep_x1_n = cep_x1;
cep_x1_n(locs(1)-1:locs(1)+1) = 0;
cep_x1_n(locs(2)-1:locs(2)+1) = 0;


logmag_X1h_dash = real(fft(cep_x1_h, N));
logmag_X1n_dash = real(fft(cep_x1_n, N));



figure(2);
subplot(311);
plot(res_voiced);
subplot(312);
plot(cep_x1);
subplot(313);
plot(res_voiced_dash);


figure(3);
subplot(511);
plot(cep_x(2:end));
subtitle("Cepstrum of voiced frame");
subplot(512);
plot(cep_x1_h);
subtitle("Harmonic part of cepstrum");
subplot(513);
plot(cep_x1_n);
subtitle("Noise part of cepstrum");
subplot(514);
plot(logmag_X1h_dash);
subtitle("Log Magnitude spectrum of harmonic part");
subplot(515);
plot(logmag_X1n_dash);
subtitle("Log Magnitude spectrum of noise part");

% Frequency bins (FFT coefficients) in noise region 'Fr'
noise_comps = zeros(length(logmag_X1n_dash),1);
noise_comps(logmag_X1n_dash >= logmag_X1h_dash) = logmag_X1n_dash(logmag_X1n_dash >= logmag_X1h_dash);

% Frequency bins (FFT coefficients) in harmonic region 'Fp'
harmonic_comps = zeros(length(logmag_X1h_dash),1);
harmonic_comps(logmag_X1n_dash < logmag_X1h_dash) = logmag_X1h_dash(logmag_X1n_dash < logmag_X1h_dash);


figure(4);
subplot(411);
plot(logmag_X1h_dash);
subtitle('Log Magnitude Spectrum: Harmonic part');
ylim([-1.8,1.8]);
subplot(412);
plot(logmag_X1n_dash);
subtitle('Log Magnitude Spectrum: Noise part');
ylim([-1.8,1.8]);
subplot(413);
plot(harmonic_comps);
subtitle('Log Magnitude Spectrum: Harmonic components');
ylim([-1.8,1.8]);
subplot(414);
plot(noise_comps);
subtitle('Log Magnitude Spectrum: Noise components');
ylim([-1.8,1.8]);

% Signal Reconstruction

mag_harmonicX_dash = exp(harmonic_comps);
res_harmonicX_dash = mag_harmonicX_dash.*exp(1i*phase_X);
res_harmonic_dash = real(ifft(res_harmonicX_dash, N));

mag_noiseX_dash = exp(noise_comps);
res_noiseX_dash = mag_noiseX_dash.*exp(1i*phase_X);
res_noise_dash = real(ifft(res_noiseX_dash, N));

figure(5);
subplot(311);
plot(res_voiced);
subtitle("Residual of voiced frame");
subplot(312);
plot(res_harmonic_dash);
subtitle("Reconstructed Periodic frame");
subplot(313);
plot(res_noise_dash);
subtitle("Reconstructed Aperiodic frame");


figure(6);
subplot(211);
cep_h = real(ifft(log(abs(fft(res_harmonic_dash, N)))));
plot(cep_h);
subplot(212);
cep_n = real(ifft(log(abs(fft(res_noise_dash, N)))));
plot(cep_n);

