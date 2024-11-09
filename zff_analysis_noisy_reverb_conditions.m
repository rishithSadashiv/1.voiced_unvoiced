% To analyse impact of voiced-unvoiced decomposition on voiced components.


% free field model: x(n) = as(n+t) + v(n)  [a-attenuation constant,
% v(n)-noise, t-tdoa]


% Room reverberent model: x(n) = h(n)*s(n) + v(n) [h(n)-channel impulse
% response]

% In noisy conditions, since noise is an additive term, it can be
% considered the ZFF impulses are robust to noise upto certain SNR


clc;
clear all;
close all;

[d,fs]=audioread('./arctic_b0399.wav');

d=d-mean(d);
d=d./max(abs(d));
original_aud = d;


[d_reverb,fs]=audioread('1.audio_withReverb_t60_400.wav');
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

[m1_epochlocs,m1_zsp1,m1_vgclocssp1] = EpochsbyZFF(original_aud, fs);
m1_imp_train = zeros(1,length(original_aud));
m1_imp_train(m1_epochlocs) = 1;
m1_imp_train = m1_imp_train(:);


figure(2);
subplot(311);
plot(original_aud);
subplot(312);
plot(m1_zsp1);
subplot(313);
plot(m1_vgclocssp1);

ms30 = 0.03*fs;
voiced_index = 5520;
unvoiced_index = 22336;
figure(3);
subplot(311);
plot(original_aud(voiced_index:voiced_index+ms30-1));
subplot(312);
plot(abs(m1_zsp1(voiced_index:voiced_index+ms30-1)));
% subplot(313);
% plot(m1_vgclocssp1(voiced_index:voiced_index+ms30-1));
a = abs(m1_zsp1(voiced_index:voiced_index+ms30-1));
disp(sum(buffer(a.^2, length(a))));


figure(4);
subplot(311);
plot(original_aud(unvoiced_index:unvoiced_index+ms30-1));
subplot(312);
plot(abs(m1_zsp1(unvoiced_index:unvoiced_index+ms30-1)));
% subplot(313);
% plot(m1_vgclocssp1(voiced_index:voiced_index+ms30-1));

b = abs(m1_zsp1(unvoiced_index:unvoiced_index+ms30-1));
disp(sum(buffer(b.^2, length(b))));




[m2_epochlocs,m2_zsp1,m2_vgclocssp1] = EpochsbyZFF(noise_aud, fs);
m2_imp_train = zeros(1,length(noise_aud));
m2_imp_train(m2_epochlocs) = 1;
m2_imp_train = m2_imp_train(:);

[m3_epochlocs,m3_zsp1,m3_vgclocssp1] = EpochsbyZFF(reverb_aud, fs);
m3_imp_train = zeros(1,length(reverb_aud));
m3_imp_train(m3_epochlocs) = 1;
m3_imp_train = m3_imp_train(:);



original_voiced = original_aud(voiced_index:voiced_index+ms30-1);
m1_imp_train_voiced = m1_imp_train(voiced_index:voiced_index+ms30-1);
noise_voiced = noise_aud(voiced_index:voiced_index+ms30-1);
m2_imp_train_voiced = m2_imp_train(voiced_index:voiced_index+ms30-1);
reverb_voiced = reverb_aud(voiced_index:voiced_index+ms30-1);
m3_imp_train_voiced = m3_imp_train(voiced_index:voiced_index+ms30-1);



figure(1);
subplot(321);
plot(original_voiced);
subplot(322);
plot(m1_imp_train_voiced);
subplot(323);
plot(noise_voiced);
subplot(324);
plot(m2_imp_train_voiced);
subplot(325);
plot(reverb_voiced);
subplot(326);
plot(m3_imp_train_voiced);




