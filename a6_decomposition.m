% Decomposition of speech signal into periodic and aperiodic parts.
% Input: Speech signal
% Output: Periodic signal, Aperiodic signal

% Based on a5_decomposition.m file


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
[res_d,lpcs]=LPres(d,fs,framesize,frameshift,lporder,1);

%% Preparation for Short-term windowing
framesize_samples = framesize/1000*fs;
frameshift_samples = frameshift/1000*fs;
nx = length(res_d);
nframes = fix((nx-framesize_samples)/frameshift_samples+1);
frames= zeros(framesize_samples,nframes);

frameindeces = 1 + (0:(nframes-1))*frameshift_samples;
rowindex = (1:framesize_samples)';
frames(:) = res_d(rowindex(:,ones(1,nframes))+frameindeces(ones(framesize_samples,1),:)-1);
win = hamming(framesize_samples);
frames_windowed = win(:,ones(1,nframes)).*frames;

%% Apply Periodic-Aperiodic decomposition

N=framesize_samples*2;

periodic_frames = zeros(framesize_samples,nframes);
aperiodic_frames = zeros(framesize_samples,nframes);

for i=1:nframes
    res_framew = frames_windowed(:,i);
    [res_harmonic_i, res_noise_i] = p_ap_decomp(res_framew,N,fs);
    periodic_frames(:,i)=res_harmonic_i;
    aperiodic_frames(:,i)=res_noise_i;
end

%% Overlap-add synthesis
periodic_lpr = overlap_add(periodic_frames, framesize_samples, frameshift_samples, nframes);
aperiodic_lpr = overlap_add(aperiodic_frames, framesize_samples, frameshift_samples, nframes);


%% Linear Prediction synthesis
LPCs = lpcs(1:nframes,:); % removing LP coefficients of the last frame
plotflag=0;
periodic_component = SynthSpeech_v5(periodic_lpr,LPCs,lporder,framesize_samples,frameshift_samples,plotflag);
aperiodic_component = SynthSpeech_v5(aperiodic_lpr,LPCs,lporder,framesize_samples,frameshift_samples,plotflag);

t=1:length(aperiodic_component);
t=t/fs;
figure();
subplot(321);
plot(t,res_d(1:length(aperiodic_lpr)));
subtitle("LPR of original signal");
subplot(323);
plot(t,aperiodic_lpr);
subtitle("Aperiodic part of LPR");
subplot(325);
plot(t,periodic_lpr);
subtitle("Periodic part of LPR");
subplot(322);
plot(t,d(1:length(periodic_component)));
subtitle("Original signal");
subplot(324);
plot(t,aperiodic_component);
subtitle("aperiodic part of original signal");
subplot(326);
plot(t,periodic_component);
subtitle("Periodic part of original signal");


%% Functions

function signal = overlap_add(frames, framesize_samples, frameshift_samples, nframes)
    signal=zeros((nframes-1)*frameshift_samples+framesize_samples,1);
    idx=(1:framesize_samples)';
    start=0;
    for l=1:nframes
       signal(start+idx)=signal(start+idx)+frames(:,l);
       start=start+frameshift_samples;
    end


end


function [res_harmonic_i, res_noise_i] = p_ap_decomp(res_framew, N,fs)
    
    % Cepstral analysis (Obtain the initial estimate of aperiodic component)
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
    res_noiseX_dash = zeros(N,1);
    res_noiseX_dash(noise_comps==1) = res_X(noise_comps==1);
    res_noise_dash = real(ifft(res_noiseX_dash, N));

    % Iterative algorithm for bandlimited signal extrapolation (to reconstruct the aperiodic part)
    % Finite duration constraint in time-domain
    % Known noise samples constraint in freq-domain
    
    % First iteration 
    R0 = res_noiseX_dash;
    r0 = ifft(R0,N);
    r0_hat = r0;
    r0_hat(N/2:N) = 0;
    rm_minus1_hat = r0_hat;

    % mth iteration
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

    % Subtracting noise estimate from original LP residual to get the harmonic estimate.
    res_noise_i = gm(1:N/2);
    res_harmonic_i = res_framew - res_noise_i;

end






