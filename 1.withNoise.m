% Computing TDoA when AWGN is present.

clc;
clear all;
close all;

[d,fs]=audioread('./arctic_b0399.wav');
d=d(:,1);

d=d-mean(d);
d=d./max(abs(d));


ms30 = 0.03 * fs;
shift = 16;
voiced_index = 11648;
unvoiced_index = 22336;
SNR = 20;

mic1_aud = 0.9*d;
mic2_aud = 0.8*d;
mic1_aud_noisy = awgn(mic1_aud, SNR, "measured");
mic2_aud_noisy = awgn(mic2_aud, SNR, "measured");

mic1_voiced = mic1_aud_noisy(voiced_index:voiced_index+ms30-1);
mic2_voiced = mic2_aud_noisy(voiced_index+shift: voiced_index+shift+ms30-1);
mic1_unvoiced = mic1_aud_noisy(unvoiced_index:unvoiced_index+ms30-1);
mic2_unvoiced = mic2_aud_noisy(unvoiced_index+shift: unvoiced_index+shift+ms30-1);

%%% Excitation Source Method

expected_delay_voiced = tdoa_excitationSource_method(mic1_aud_noisy, mic2_aud_noisy, voiced_index, ms30, fs, shift);
expected_delay_unvoiced = tdoa_excitationSource_method(mic1_aud_noisy, mic2_aud_noisy, unvoiced_index, ms30, fs, shift);

disp("Excitation Source Method:")
disp("TDoA (voiced): " + expected_delay_voiced);
disp("TDoA (unvoiced): " + expected_delay_unvoiced);


%%% GCC-PHAT method
[estimated_delay_voiced, cc1, peak_position1, peak_value1] = delay_gcc(mic1_voiced, mic2_voiced, "cc");
[estimated_delay_unvoiced, cc2, peak_position2, peak_value2] = delay_gcc(mic1_unvoiced, mic2_unvoiced, "cc");
disp("GCC-PHAT Method:")
disp("TDoA (voiced): " + expected_delay_voiced);
disp("TDoA (unvoiced): " + expected_delay_unvoiced);










