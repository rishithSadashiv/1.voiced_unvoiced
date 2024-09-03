%%


clc;
clear all;
close all;

[wav1,fs]=audioread('./20240901/15cm/30degree/mic1_15cm_30deg.WAV');
[wav2,fs]=audioread('./20240901/15cm/30degree/mic2_15cm_30deg.WAV');

% wav1=resample(wav1,16000,fs);
% wav2=resample(wav2,16000,fs);
% fs=16000;



wav1=wav1(:,2);
wav2=wav2(:,4);
%d=resample(d,8000,fs);
%fs=16000;
wav1=wav1-mean(wav1);
wav1=wav1./max(abs(wav1));

wav2 = wav2-mean(wav2);
wav2=wav2./max(abs(wav2));


ms30 = 0.03 * fs;
voiced_index = 25.00*fs;
SNR = 10;


mic1_voiced = wav1(voiced_index:voiced_index+ms30-1);
mic2_voiced = wav2(voiced_index:voiced_index+ms30-1);

hamming_w = hamming(length(mic1_voiced));
mic1_voiced_w = mic1_voiced.*hamming_w;
mic2_voiced_w = mic2_voiced.*hamming_w;

mic1_aud_lpr = LPres(wav1,fs,20,5,10,1);
mic2_aud_lpr = LPres(wav2,fs,20,5,10,1);
% 
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Using Excitation Source Information
% 
mic1_voiced_lpr = mic1_aud_lpr(voiced_index:voiced_index+ms30-1);
mic2_voiced_lpr = mic2_aud_lpr(voiced_index: voiced_index+ms30-1);
% mic1_voiced_lpr = mic1_voiced_lpr.*hamming_w;
% mic2_voiced_lpr = mic2_voiced_lpr.*hamming_w;
% 
mic1_voiced_lpr_henv = HilbertEnv(mic1_voiced_lpr,fs);
mic2_voiced_lpr_henv = HilbertEnv(mic2_voiced_lpr, fs);
% 
% 
figure(1);
a(1) = subplot(321);
plot(mic1_voiced,'k');
title('Mic1 Voiced');
a(2) = subplot(322);
plot(mic2_voiced,'k');
title('Mic2 Voiced');

a(3) = subplot(323);
plot(mic1_voiced_lpr);
title('Mic1 LPR');

a(4) = subplot(324);
plot(mic2_voiced_lpr);
title('Mic2 LPR');

a(5) = subplot(325);
plot(mic1_voiced_lpr_henv);
title('Mic1 LPR Hilbert Env');

a(6) = subplot(326);
plot(mic2_voiced_lpr_henv);
title('Mic2 LPR Hilbert Env');


% Z1 = fft(mic1_voiced_lpr_henv,ms30);
% Z2 = fft(mic2_voiced_lpr_henv,ms30);
% % 
% % figure(2);
% % a(1) = subplot(2,2,1);
% % plot(abs(Z1));
% % a(2) = subplot(2,2,2);
% % plot(angle(Z1));
% % a(3) = subplot(2,2,3);
% % plot(abs(Z2));
% % a(4) = subplot(2,2,4);
% % plot(angle(Z2))
% 
disp("Using Excitation Information");
[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_voiced_lpr_henv, mic2_voiced_lpr_henv, "cc");
disp("Voiced segment: estimated delay = " + estimated_delay)
% 
figure(2);
plot(cc)
title('H Env CCF of voiced segment');
hold on;
scatter(peak_position, peak_value, "filled");
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
unvoiced_index = 24.320*fs;
mic1_unvoiced = wav1(unvoiced_index:unvoiced_index+ms30-1);
mic2_unvoiced = wav2(unvoiced_index:unvoiced_index+ms30-1);
% mic1_unvoiced = mic1_unvoiced.*hamming_w;
% mic2_unvoiced = mic2_unvoiced.*hamming_w;
% 
mic1_unvoiced_lpr = mic1_aud_lpr(unvoiced_index:unvoiced_index+ms30-1);
mic2_unvoiced_lpr = mic2_aud_lpr(unvoiced_index: unvoiced_index+ms30-1);
% mic1_unvoiced_lpr = mic1_unvoiced_lpr.*hamming_w;
% mic2_unvoiced_lpr = mic2_unvoiced_lpr.*hamming_w;
mic1_unvoiced_lpr_henv = HilbertEnv(mic1_unvoiced_lpr, fs);
mic2_unvoiced_lpr_henv = HilbertEnv(mic2_unvoiced_lpr, fs);

figure(3);
a(1) = subplot(321);
plot(mic1_unvoiced,'k');
title('Mic1 Unvoiced');
a(2) = subplot(322);
plot(mic2_unvoiced,'k');
title('Mic2 Unvoiced');

a(3) = subplot(323);
plot(mic1_unvoiced_lpr);
title('Mic1 LPR');
a(4) = subplot(324);
plot(mic2_unvoiced_lpr);
title('Mic2 LPR')

a(5) = subplot(325);
plot(mic1_unvoiced_lpr_henv);
title('Mic1 LPR Hilbert Env')

a(6) = subplot(326);
plot(mic2_unvoiced_lpr_henv);
title('Mic2 LPR Hilbert Env')

[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_unvoiced_lpr_henv, mic2_unvoiced_lpr_henv, "cc");
disp("Unvoiced segment: estimated delay = " + estimated_delay)

figure(4);
plot(cc);
title('H Env CCF of Unvoiced ')
hold on;
scatter(peak_position, peak_value, "filled");
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%% Using GCC-PHAT
% 
disp("Using GCC-PHAT");
[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_voiced_w, mic2_voiced_w, "phat");
disp("Voiced segment: estimated delay = " + estimated_delay)
% 
% 
figure(5);
plot(cc)
title('CCF GCC-PHAT voiced')

hold on;
scatter(peak_position, peak_value, "filled");
% 
[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_unvoiced, mic2_unvoiced, "phat");
disp("Unvoiced segment: estimated delay = " + estimated_delay)

figure(6);
plot(cc)
title('CCF GCC-PHAT unvoiced')
hold on;
scatter(peak_position, peak_value, "filled");
% 
%

figure(7);
a(1) = subplot(211);
plot(mic1_voiced,'k');
a(2) = subplot(212);
plot(mic2_voiced,'k');

figure(8);
a(1) = subplot(211);
plot(mic1_unvoiced,'k');
a(2) = subplot(212);
plot(mic2_unvoiced,'k');
