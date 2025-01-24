% Decomposition of speech signal into periodic and aperiodic parts.
% Input: Speech frame
% Output: Periodic part of LP Residual frame, Aperiodic part of LP Residual frame

% Based on a2_decomposition.m file


clc;
clear all;
close all;

[d,fs]=audioread('./arctic_b0399.wav');
d=d-mean(d);
d=d./max(abs(d));

%% LP residual 
% parameters
framesize = 32;
frameshift = 4;
lporder = 10;
res_d=LPres(d,fs,framesize,frameshift,lporder,1);

% Working with one frame
voiced_index = 5520;%11648;
ms32 = 0.032*fs;
res_frame = res_d(voiced_index:voiced_index+ms32-1);
res_framew = res_frame.*hamming(length(res_frame));

%% Cepstral analysis (Obtain the initial estimate of aperiodic component)
N = 1024;
res_X = fft(res_framew, N);
mag_X = abs(res_X);
phase_X = angle(res_X);
logmag_X = log(mag_X);
cep_x = ifft(logmag_X, N);

% Removing vocal tract system characteristics.
cep_x1 = cep_x;
ltf_sample = 0.002*fs; % Low time lifter cutoff 0.002 s
cep_x1(1:ltf_sample) = 0;
cep_x1(end-ltf_sample+1:end) = 0;

% Keeping the harmonic part
[vals, locs] = maxk(cep_x1, 2);
cep_x1_h = zeros(length(cep_x1),1);
cep_x1_h(locs(1)-2:locs(1)+2) = cep_x1(locs(1)-2:locs(1)+2);
cep_x1_h(locs(2)-2:locs(1)+2) = cep_x1(locs(2)-2:locs(1)+2);

% Removing the harmonic part
cep_x1_n = cep_x1;
cep_x1_n(locs(1)-2:locs(1)+2) = 0;
cep_x1_n(locs(2)-2:locs(2)+2) = 0;

% Applying FFT to the liftered signal
mag_X1h_dash = exp(real(fft(cep_x1_h, N)));
mag_X1n_dash = exp(real(fft(cep_x1_n, N)));

% Frequency bins (FFT coefficients) in noise region 'Fr'
noise_comps = zeros(length(mag_X1n_dash),1);
noise_comps(mag_X1n_dash >= mag_X1h_dash) = 1; % indicating the indeces which correspond to noise components


% Signal Reconstruction

% res_noiseX_dash = res_X .* noise_comps;
res_noiseX_dash = zeros(N,1);
res_noiseX_dash(noise_comps==1) = res_X(noise_comps==1);
res_noise_dash = real(ifft(res_noiseX_dash, N));


%% Iterative algorithm for bandlimited signal extrapolation (to reconstruct the aperiodic part)
% Finite duration constraint in time-domain
% known noise samples constraint in freq-domain

% First iteration 

R0 = res_noiseX_dash;
r0 = ifft(R0,N);
r0_hat = r0;
r0_hat(N/2:N) = 0;

rm_minus1_hat = r0_hat;

%% mth iteration

i=1;
n_iterations = 10;
while (i <= n_iterations)

Rm_minus1_hat = fft(rm_minus1_hat,N);
Rm = Rm_minus1_hat;
Rm(noise_comps==1) = res_X(noise_comps==1);
gm = ifft(Rm,N);
rm_minus1_hat = gm;
rm_minus1_hat(N/2:N)=0;
i = i+1;

end

%% Subtracting noise estimate from original LP residual to get the harmonic estimate.

res_noise_i = gm(1:N/2);
res_harmonic_i = res_framew - res_noise_i;

figure();
subplot(321);
plot(res_framew);
title("Waveform plots");
subtitle("Windowed LP residual");
subplot(323);
plot(res_noise_i);
subtitle("Aperiodic component of LPR (iterative method estimate)");
subplot(325);
plot(res_harmonic_i);
subtitle("Periodic component of LPR (iterative method estimate)");
subplot(322);
plot(log(abs(res_X(1:N/2))));
subtitle("Windowed LPR frame");
title("Log magnitude spectrum")
subplot(324);
a1=log(abs(fft(res_noise_i,N)));
plot(a1(1:N/2));
subtitle("Aperiodic component of LPR");
subplot(326);
a2=log(abs(fft(res_harmonic_i,N)));
plot(a2(1:N/2));
subtitle("Aperiodic component of LPR");


figure();
subplot(222);
plot(res_noise_dash(1:N/2));
subtitle("Cepstral analysis estimate");
title("Noise Component of LPR: Waveform");
subplot(224);
plot(res_noise_i);
subtitle("Iterative method estimate");
subplot(221);
a0 = log(abs(res_noiseX_dash+0.00001));  % mag_spec has zeros. log zero becomes -inf. To avoid this, small number is added.
plot(a0(1:N/2));
subtitle("Cepstral analysis estimate");
title("Noise component of LPR: Log Magnitude Spectrum");
subplot(223);
plot(a1(1:N/2));
subtitle("Iterative methods esimtate");




