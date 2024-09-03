% gcc_phat other implementation

clc;
clear all;
close all;

% [d,fs]=audioread('./arctic_b0399.wav');
% d=d(:,1);
% 
% d=d-mean(d);
% d=d./max(abs(d));

[wav1,fs]=audioread('./20240901/15cm/30degree/mic1_15cm_30deg.WAV');
[wav2,fs]=audioread('./20240901/15cm/30degree/mic2_15cm_30deg.WAV');
% start_index = 25.141*fs;
% end_index = 24.320*fs;

wav1=wav1(:,2);
wav2=wav2(:,4);
%d=resample(d,8000,fs);
%fs=16000;
wav1=wav1-mean(wav1);
wav1=wav1./max(abs(wav1));

wav2 = wav2-mean(wav2);
wav2=wav2./max(abs(wav2));



ms30 = 0.03 * fs;
shift = 16;
% voiced_index = 11648;
voiced_index = 25.141*fs;
% unvoiced_index = 22336;
unvoiced_index = 24.320*fs;
% voiced_index = 16000;

% mic1_aud = 0.9*d;
% mic2_aud = 0.8*d;

mic1_voiced = wav1(voiced_index:voiced_index+ms30-1);
mic2_voiced = wav2(voiced_index:voiced_index+ms30-1);
wlen = length(mic1_voiced);
w = hann(wlen);
mic1_voiced = mic1_voiced.*w;
mic2_voiced = mic2_voiced.*w;

Nfft = wlen;
GCC = fftshift(real(ifft(exp(1i*angle(fft(mic1_voiced,Nfft).*conj(fft(mic2_voiced,Nfft))))))); 

GCCn = GCC/max(GCC);

[value, peak_positions] = max(GCC);
peak_position = max(peak_positions);
peak_value = value;
estimated_delay = peak_position - length(mic1_voiced);
disp(estimated_delay);

figure(1);
plot(GCCn)
