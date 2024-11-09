% 


clc;
clear all;
close all;

[d,fs]=audioread('./arctic_b0399.wav');
% [d,fs]=audioread('new_wavs/arctic_a0023.wav');
d=d(:,1);


% For frame-wise operation
n=300;
frameshift = 0.02*fs;
shift=16;
ms30=0.03*fs;
audiolen = length(d);

zcr_array = [];
ste_array = [];

SNRs = [50,20,10,5,0,-5];
ste_array = csvread('csv_forPlots\ste_array_noNoise_2.csv');
ste_array=ste_array(:);

d_v = zeros(length(SNRs), 160);
d_uv = zeros(length(SNRs), 160);

delay_array = zeros(length(SNRs), 160);
for p=1:length(SNRs)
    
    SNR = SNRs(p);
    % disp(SNR);
    if(SNR == 50)
        d=d-mean(d);
        d=d./max(abs(d));
        mic1_aud = 0.9*d;
        mic2_aud = 0.8*d;
    else
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
    end

    [m1_epochlocs,m1_zsp1,m1_vgclocssp1] = EpochsbyZFF(mic1_aud, fs);
    m1_imp_train = zeros(1,length(mic1_aud));
    m1_imp_train(m1_epochlocs) = 1;
    m1_imp_train = m1_imp_train(:);
    
    [m2_epochlocs,m2_zsp1,m2_vgclocssp1] = EpochsbyZFF(mic2_aud, fs);
    m2_imp_train = zeros(1,length(mic2_aud));
    m2_imp_train(m2_epochlocs) = 1;
    m2_imp_train = m2_imp_train(:);

    computed_delays = zeros(n,1);
    v=zeros(n,1);
    u=zeros(n,1);
    for i = 1:n
        start1 = ((i-1)*frameshift)+1;
        start2 = start1+shift;
        end1 = start1+ms30-1;
        end2 = end1+shift;
        if(end2>audiolen)
            break;
        end
        
        mic1_seg = mic1_aud(start1:end1);
        % zcr_frame = zerocrossrate(mic1_seg);
        ste_frame1 = sum(buffer(mic1_seg.^2, length(mic1_seg)));
    
        % zcr_array(i) = zcr_frame;
        ste_array1(i) = ste_frame1;
        ste_frame = ste_array(i);

        if(ste_frame > 0.15)

            mic1_seg = m1_imp_train(start1:end1);
            mic2_seg = m2_imp_train(start2:end2);
            
            estimated_delay = computeDelayManual(mic1_seg, mic2_seg);
            computed_delays(i) = estimated_delay;
        
        else
            mic1_seg = mic1_aud(start1:end1);
            mic2_seg = mic2_aud(start2:end2);
            estimated_delay = gcc_phat_alternative(mic1_seg, mic2_seg);
            computed_delays(i) = estimated_delay;
        end
        ste_frame = ste_array(i);
        if(ste_frame > 0.15)
            v(i) = estimated_delay;
        else
            u(i)=estimated_delay;
        end
    end
    delay_array(p,1:160) = computed_delays(1:160);
    d_v(p,:)=v(1:160);
    d_uv(p,:)=u(1:160);
end

writematrix(delay_array, 'csv_forPlots/merged_delay_array_arctic_b0399.csv')
writematrix(d_v, 'csv_forPlots/merged_delay_vowel_arctic_b0399.csv')

frames_index = 1:160;

figure(1);
subplot(7,1,1);
scatter(frames_index,delay_array(1,:), 16,'filled');
ylim([-50 50])
subplot(712);
scatter(frames_index, delay_array(2,:), 16, 'filled');
ylim([-50 50])
subplot(713);
scatter(frames_index, delay_array(3,:), 16, 'filled');
ylim([-50 50])
subplot(714);
scatter(frames_index, delay_array(4,:), 16, 'filled');
ylim([-50 50])
subplot(715);
scatter(frames_index, delay_array(5,:), 16,'filled');
ylim([-50 50])
subplot(716);
scatter(frames_index, delay_array(6,:), 16,'filled');
ylim([-50 50])
subplot(717);
plot(ste_array(1:160));

% 
% avg_delay_20dB = mean(delay_array(2,:));
% avg_delay_5dB = mean(delay_array(4,:));
% avg_delay_neg5dB = mean(delay_array(6,:));
% 
% disp("Avg 20dB:"+avg_delay_20dB);
% disp("Avg 5dB:"+avg_delay_5dB);
% disp("Avg -5dB:"+avg_delay_neg5dB);
% 
% correctCount_20dB = sum(delay_array(2,:)==16);
% correctCount5dB = sum(delay_array(4,:)==16);
% correctCoungNeg5dB = sum(delay_array(6,:)==16);
% 
% disp("Correct count 20dB:"+correctCount_20dB);
% disp("Correct count 5dB:"+correctCount5dB);
% disp("Correct count -5dB:"+correctCoungNeg5dB);
% 
% var_20dB = var(delay_array(2,:));
% var_5dB = var(delay_array(4,:));
% var_Neg5dB = var(delay_array(6,:));
% 
% disp("Variance 20dB:"+ var_20dB);
% disp("Variance 5dB:"+ var_5dB);
% disp("Variance -5dB:"+ var_Neg5dB);



correctCount_20dB = sum(d_v(2,:)==16);
correctCount5dB = sum(d_v(4,:)==16);
correctCoungNeg5dB = sum(d_v(6,:)==16);

disp("Correct count 20dB:"+correctCount_20dB);
disp("Correct count 5dB:"+correctCount5dB);
disp("Correct count -5dB:"+correctCoungNeg5dB);

correctCount_20dB = sum(d_uv(2,:)==16);
correctCount5dB = sum(d_uv(4,:)==16);
correctCoungNeg5dB = sum(d_uv(6,:)==16);

disp("Correct count 20dB:"+correctCount_20dB);
disp("Correct count 5dB:"+correctCount5dB);
disp("Correct count -5dB:"+correctCoungNeg5dB);

nbins = 100;
figure(2);
% boxchart(delay_array');
subplot(611);
histogram(delay_array(1,:),nbins);
subtitle('No noise');
xlim([-100 100])
subplot(612);
histogram(delay_array(2,:),nbins);
subtitle('20dB');
xlim([-100 100])
subplot(613);
histogram(delay_array(3,:),nbins);
subtitle('10dB')
xlim([-100 100])
subplot(614);
histogram(delay_array(4,:),nbins);
subtitle('5dB');
xlim([-100 100])
subplot(615);
histogram(delay_array(5,:),nbins);
subtitle('0dB');
xlim([-100 100])
subplot(616);
histogram(delay_array(6,:),nbins);
subtitle('-5dB')
xlim([-100 100])




function delay = computeDelayManual(m1_imp, m2_imp)
    
    m1=find(m1_imp);
    m2=find(m2_imp);
    p = length(m1);
    q = length(m2);
    if(p~=q)
        if(p>q)
            i1=m1;
            i2=m2;
        elseif(p<q)
            i1=m2;
            i2=m1;
        end
        a=abs(i1(1)-i2(1));
        b=abs(i1(2)-i2(1));
        if(a>b) % Incase of additional impulse appearing
            i1 = i1(2:end);
        else
            i1 = i1(1:end-1);
        end
        if(length(i1) ~= length(i2)) % Incase of unvoiced segments (random impulses)
            i1 = i1(1:end-(length(i1)- length(i2)));
        end
        if(p<q)
            x=i1;
            i1=i2;
            i2=x;
        end
    
    else
        i1=m1;
        i2=m2;
    end
    d=[];

    for s = 1:length(i1) 
        d(s)=i1(s)-i2(s);
    end
    delay=round(mean(d));

end


