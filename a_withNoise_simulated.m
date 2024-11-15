% Computing TDoA when AWGN is present.

clc;
clear all;
close all;

[d,fs]=audioread('./arctic_b0399.wav');
% [d,fs]=audioread('./car_sn10.wav');
d=d(:,1);

d=d-mean(d);
d=d./max(abs(d));


ms30 = 0.03 * fs;
shift = 16;
% voiced_index = 11648;
unvoiced_index = 22336;
voiced_index = 16000;


SNR = 30;

mic1_aud = 0.9*d;
mic2_aud = 0.8*d;

% Add noise: method 1
% mic1_aud_noisy = awgn(mic1_aud, SNR, "measured");
% mic2_aud_noisy = awgn(mic2_aud, SNR, "measured");

% Add noise: method 2
Px = 0.5*(mean([mic1_aud(:); mic2_aud(:)].^2));
Pn = Px*10^(-SNR/10);

rng('default')
mic1_aud_noisy = mic1_aud + sqrt(Pn)*randn(size(mic1_aud));
mic2_aud_noisy = mic2_aud + sqrt(Pn)*randn(size(mic2_aud));


mic1_voiced = mic1_aud_noisy(voiced_index:voiced_index+ms30-1);
mic2_voiced = mic2_aud_noisy(voiced_index+shift: voiced_index+shift+ms30-1);
hamming_w = hamming(length(mic1_voiced));
mic1_voiced = mic1_voiced.*hamming_w;
mic2_voiced = mic2_voiced.*hamming_w;

mic1_unvoiced = mic1_aud_noisy(unvoiced_index:unvoiced_index+ms30-1);
mic2_unvoiced = mic2_aud_noisy(unvoiced_index+shift: unvoiced_index+shift+ms30-1);
mic1_unvoiced = mic1_unvoiced.*hamming_w;
mic2_unvoiced = mic2_unvoiced.*hamming_w;

%%% Excitation Source Method

estimated_delay_voiced1 = tdoa_excitationSource_method(mic1_aud_noisy, mic2_aud_noisy, voiced_index, ms30, fs, shift);
estimated_delay_unvoiced1 = tdoa_excitationSource_method(mic1_aud_noisy, mic2_aud_noisy, unvoiced_index, ms30, fs, shift);

disp("Excitation Source Method:")
disp("TDoA (voiced): " + samples_to_time(estimated_delay_voiced1,fs));
disp("TDoA (unvoiced): " + samples_to_time(estimated_delay_unvoiced1,fs));


%%% GCC-PHAT method
[estimated_delay_voiced2, cc1, peak_position1, peak_value1] = delay_gcc(mic1_voiced, mic2_voiced, "phat");
[estimated_delay_unvoiced2, cc2, peak_position2, peak_value2] = delay_gcc(mic1_unvoiced, mic2_unvoiced, "phat");
disp("GCC-PHAT Method:")
disp("TDoA (voiced): " + samples_to_time(estimated_delay_voiced2, fs));
disp("TDoA (unvoiced): " + samples_to_time(estimated_delay_unvoiced2, fs));










