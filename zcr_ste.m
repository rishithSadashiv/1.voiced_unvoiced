%



clc;
clear all;
close all;

% [d,fs]=audioread('./arctic_b0399.wav');
[d,fs]=audioread('new_wavs/arctic_a0023.wav');
d=d(:,1);
noise = 0;

if(noise == 1)
    
    SNR = -5;
    Px = 0.5*(mean([d(:)].^2));
    Pn = Px*10^(-SNR/10);
    
    rng('default');
    d1 = d + sqrt(Pn)*randn(size(d));
    d2 = d + sqrt(Pn)*randn(size(d));
    
    d1=d1-mean(d1);
    d1=d1./max(abs(d1));
    
    d2=d2-mean(d2);
    d2=d2./max(abs(d2));
    
    mic1_aud = 0.9*d1;
    mic2_aud = 0.8*d2;


else
    d=d-mean(d);
    d=d./max(abs(d));
    mic1_aud = 0.9*d;
    mic2_aud = 0.8*d;
end


n=300;
frameshift = 0.02*fs;
shift=16;
ms30=0.03*fs;
audiolen = length(d);

zcr_array = [];
ste_array = [];
for i = 1:n
    start1 = ((i-1)*frameshift)+1;
    start2 = start1+shift;
    end1 = start1+ms30-1;
    end2 = end1+shift;
    if(end2>audiolen)
        break;
    end
    mic1_seg = mic1_aud(start1:end1);
    % mic2_seg = mic2_aud(start2:end2);
    
    zcr_array(i) = zerocrossrate(mic1_seg);
    ste_array(i) = sum(buffer(mic1_seg.^2, length(mic1_seg)));

end

disp(max(ste_array));

ste_array = ste_array./max(ste_array);

figure(1);
subplot(311);
plot(d);
subplot(312);
plot(zcr_array(:))
subplot(313);
plot(ste_array(:))

writematrix(ste_array, 'csv_forPlots/ste_array_noNoise_2.csv');
    
