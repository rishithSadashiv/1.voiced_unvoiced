
clc;
clear all;
close all;

% [d,fs]=audioread('./arctic_a0001.wav');
[d,fs]=audioread('./arctic_b0399.wav');
d=d(:,1);
%d=resample(d,8000,fs);
%fs=16000;
d=d-mean(d);
d=d./max(abs(d));

ms30 = 0.03 * fs;
shift = 16;
% voiced_index = 4500;
voiced_index = 11648;
% voiced = d(4400:4400+2*ms30);

mic1_aud = 0.9*d;
mic2_aud = 0.8*d;

mic1_voiced = mic1_aud(voiced_index:voiced_index+ms30-1);
mic2_voiced = mic2_aud(voiced_index+shift: voiced_index+shift+ms30-1);

mic1_aud_lpr = LPres(mic1_aud,fs,20,5,10,1);
mic2_aud_lpr = LPres(mic2_aud,fs,20,5,10,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Using Excitation Source Information

mic1_voiced_lpr = mic1_aud_lpr(voiced_index:voiced_index+ms30-1);
mic2_voiced_lpr = mic2_aud_lpr(voiced_index+shift: voiced_index+shift+ms30-1);
mic1_voiced_lpr_henv = HilbertEnv(mic1_voiced_lpr,fs);
mic2_voiced_lpr_henv = HilbertEnv(mic2_voiced_lpr, fs);


figure(1);
a(1) = subplot(321);
plot(mic1_voiced,'k');
a(2) = subplot(322);
plot(mic2_voiced,'k');

a(3) = subplot(323);
plot(mic1_voiced_lpr);
a(4) = subplot(324);
plot(mic2_voiced_lpr);

a(5) = subplot(325);
plot(mic1_voiced_lpr_henv);
a(6) = subplot(326);
plot(mic2_voiced_lpr_henv);


% Z1 = fft(mic1_voiced_lpr_henv,ms30);
% Z2 = fft(mic2_voiced_lpr_henv,ms30);
% 
% figure(2);
% a(1) = subplot(2,2,1);
% plot(abs(Z1));
% a(2) = subplot(2,2,2);
% plot(angle(Z1));
% a(3) = subplot(2,2,3);
% plot(abs(Z2));
% a(4) = subplot(2,2,4);
% plot(angle(Z2))

[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_voiced_lpr_henv, mic2_voiced_lpr_henv, "cc");
disp("Voiced segment: estimated delay = " + estimated_delay)

figure(2);
plot(cc)
hold on;
scatter(peak_position, peak_value, "filled");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unvoiced_index = 44000;
unvoiced_index = 22336;
mic1_unvoiced = mic1_aud(unvoiced_index:unvoiced_index+ms30-1);
mic2_unvoiced = mic2_aud(unvoiced_index+shift: unvoiced_index+shift+ms30-1);

mic1_unvoiced_lpr = mic1_aud_lpr(unvoiced_index:unvoiced_index+ms30-1);
mic2_unvoiced_lpr = mic2_aud_lpr(unvoiced_index+shift: unvoiced_index+shift+ms30-1);
mic1_unvoiced_lpr_henv = HilbertEnv(mic1_unvoiced_lpr, fs);
mic2_unvoiced_lpr_henv = HilbertEnv(mic2_unvoiced_lpr, fs);

figure(3);
a(1) = subplot(321);
plot(mic1_unvoiced,'k');
a(2) = subplot(322);
plot(mic2_unvoiced,'k');

a(3) = subplot(323);
plot(mic1_unvoiced_lpr);
a(4) = subplot(324);
plot(mic2_unvoiced_lpr);

a(5) = subplot(325);
plot(mic1_unvoiced_lpr_henv);
a(6) = subplot(326);
plot(mic2_unvoiced_lpr_henv);


[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_unvoiced_lpr_henv, mic2_unvoiced_lpr_henv, "cc");
disp("Unvoiced segment: estimated delay = " + estimated_delay)

figure(4);
plot(cc)
hold on;
scatter(peak_position, peak_value, "filled");


%%%%%%%%%%%%%%%%%%%%%%%%%% Using GCC-PHAT

[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_voiced, mic2_voiced, "cc");
disp("Voiced segment: estimated delay = " + estimated_delay)

figure(5);
plot(cc)
hold on;
scatter(peak_position, peak_value, "filled");

[estimated_delay, cc, peak_position, peak_value] = delay_gcc(mic1_unvoiced, mic2_unvoiced, "cc");
disp("Unvoiced segment: estimated delay = " + estimated_delay)

figure(6);
plot(cc)
hold on;
scatter(peak_position, peak_value, "filled");



figure(7);
a(1) = subplot(211);
plot(mic1_unvoiced,'k');
a(2) = subplot(212);
plot(mic2_unvoiced,'k');
