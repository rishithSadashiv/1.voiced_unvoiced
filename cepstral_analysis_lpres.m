% Cepstral region in reverb condition


clc;
clear all;
close all;



[d,fs]=audioread('./arctic_b0399.wav');

d=d-mean(d);
d=d./max(abs(d));
original_aud = d;


[d_reverb,fs]=audioread('1.audio_withReverb_t60_100.wav');
d_reverb=d_reverb-mean(d_reverb);
d_reverb=d_reverb./max(abs(d_reverb));
reverb_aud = d_reverb(:,3);
reverb_aud = reverb_aud(:);


SNR = 10;
Px = 0.5*(mean([d(:)].^2));
Pn = Px*10^(-SNR/10);

rng('default');
d1 = d + sqrt(Pn)*randn(size(d));
d1=d1-mean(d1);
d1=d1./max(abs(d1));
noise_aud = d1;

%%%
ms30 = 0.03*fs;
voiced_index = 5520;
original_voiced = original_aud(voiced_index:voiced_index+ms30-1);

res_orig=LPres(original_aud,fs,20,10,10,1);
res_orig_voiced = res_orig(voiced_index:voiced_index+ms30-1);


reverb_voiced = reverb_aud(voiced_index:voiced_index+ms30-1);
res_reverb=LPres(reverb_aud,fs,20,10,10,1);
res_reverb_voiced = res_reverb(voiced_index:voiced_index+ms30-1);




figure(1);
subplot(421);
plot(original_voiced);
subplot(423);
[cep_orig_voiced,ym] = rceps(original_voiced);
plot(cep_orig_voiced(2:ms30/2));
subplot(425);
plot(res_orig_voiced);
subplot(427);
[cep_res_orig_voiced,ym] = rceps(res_orig_voiced);
% cep_res_orig_voiced=cep_res_orig_voiced./(1.01*max(abs(cep_res_orig_voiced)));
% cep_res_orig_voiced = (cep_res_orig_voiced - min(cep_res_orig_voiced))/(max(cep_res_orig_voiced) - min(cep_res_orig_voiced));
plot(cep_res_orig_voiced(2:ms30/2));

subplot(422);
plot(reverb_voiced);
subplot(424);
[cep_reverb_voiced,ym] = rceps(reverb_voiced);
plot(cep_reverb_voiced(2:ms30/2));
subplot(426);
plot(res_reverb_voiced);
subplot(428);
[cep_res_reverb_voiced,ym] = rceps(res_reverb_voiced);
plot(cep_res_reverb_voiced(2:ms30/2));



[d_reverb,fs]=audioread('1.audio_withReverb_t60_050.wav');
d_reverb=d_reverb-mean(d_reverb);
d_reverb=d_reverb./max(abs(d_reverb));
reverb_aud = d_reverb(:,3);
reverb_aud_050 = reverb_aud(:);

[d_reverb,fs]=audioread('1.audio_withReverb_t60_100.wav');
d_reverb=d_reverb-mean(d_reverb);
d_reverb=d_reverb./max(abs(d_reverb));
reverb_aud = d_reverb(:,3);
reverb_aud_100 = reverb_aud(:);

[d_reverb,fs]=audioread('1.audio_withReverb_t60_150.wav');
d_reverb=d_reverb-mean(d_reverb);
d_reverb=d_reverb./max(abs(d_reverb));
reverb_aud = d_reverb(:,3);
reverb_aud_150 = reverb_aud(:);

[d_reverb,fs]=audioread('1.audio_withReverb_t60_200.wav');
d_reverb=d_reverb-mean(d_reverb);
d_reverb=d_reverb./max(abs(d_reverb));
reverb_aud = d_reverb(:,3);
reverb_aud_200 = reverb_aud(:);

[d_reverb,fs]=audioread('1.audio_withReverb_t60_250.wav');
d_reverb=d_reverb-mean(d_reverb);
d_reverb=d_reverb./max(abs(d_reverb));
reverb_aud = d_reverb(:,3);
reverb_aud_250 = reverb_aud(:);

[d_reverb,fs]=audioread('1.audio_withReverb_t60_400.wav');
d_reverb=d_reverb-mean(d_reverb);
d_reverb=d_reverb./max(abs(d_reverb));
reverb_aud = d_reverb(:,3);
reverb_aud400 = reverb_aud(:);

reverb_voiced_050 = reverb_aud_050(voiced_index:voiced_index+ms30-1);
res_reverb_050=LPres(reverb_aud_050,fs,20,10,10,1);
res_reverb_voiced_050 = res_reverb_050(voiced_index:voiced_index+ms30-1);

reverb_voiced_100 = reverb_aud_100(voiced_index:voiced_index+ms30-1);
res_reverb_100=LPres(reverb_aud_100,fs,20,10,10,1);
res_reverb_voiced_100 = res_reverb_100(voiced_index:voiced_index+ms30-1);

reverb_voiced_150 = reverb_aud_150(voiced_index:voiced_index+ms30-1);
res_reverb_150=LPres(reverb_aud_150,fs,20,10,10,1);
res_reverb_voiced_150 = res_reverb_150(voiced_index:voiced_index+ms30-1);

reverb_voiced_200 = reverb_aud_200(voiced_index:voiced_index+ms30-1);
res_reverb_200=LPres(reverb_aud_200,fs,20,10,10,1);
res_reverb_voiced_200 = res_reverb_200(voiced_index:voiced_index+ms30-1);

reverb_voiced_250 = reverb_aud_250(voiced_index:voiced_index+ms30-1);
res_reverb_250=LPres(reverb_aud_250,fs,20,10,10,1);
res_reverb_voiced_250 = res_reverb_250(voiced_index:voiced_index+ms30-1);


[cep_orig_voiced,ym] = rceps(original_voiced);
[cep_res_reverb_voiced_050,ym] = rceps(res_reverb_voiced_050);
[cep_res_reverb_voiced_100,ym] = rceps(res_reverb_voiced_100);
[cep_res_reverb_voiced_150,ym] = rceps(res_reverb_voiced_150);
[cep_res_reverb_voiced_200,ym] = rceps(res_reverb_voiced_200);
[cep_res_reverb_voiced_250,ym] = rceps(res_reverb_voiced_250);

figure(2);
subplot(621);
plot(original_voiced);
subplot(623);
plot(reverb_voiced_050);
subplot(625);
plot(reverb_voiced_100);
subplot(627);
plot(reverb_voiced_150);
subplot(629);
plot(reverb_voiced_200);
subplot(6,2,11);
plot(reverb_voiced_250);

subplot(622);
plot(cep_orig_voiced(2:ms30/2));
subplot(624);
plot(cep_res_reverb_voiced_050(2:ms30/2));
subplot(626);
plot(cep_res_reverb_voiced_100(2:ms30/2));
subplot(628);
plot(cep_res_reverb_voiced_150(2:ms30/2));
subplot(6,2,10);
plot(cep_res_reverb_voiced_200(2:ms30/2));
subplot(6,2,12);
plot(cep_res_reverb_voiced_250(2:ms30/2));




