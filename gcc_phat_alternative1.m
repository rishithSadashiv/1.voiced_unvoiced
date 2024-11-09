% 

clc;
clear all;
close all;




[d,fs]=audioread('./arctic_b0399.wav');
% [d,fs]=audioread('new_wavs/arctic_a0023.wav');
d=d(:,1);

noise = 1;
    
if(noise == 1)
    
    SNR = 30;
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



ms30=0.03*fs;
shift=16;

voiced_index = 11648;

mic1_voiced = mic1_aud(voiced_index:voiced_index+ms30-1);
mic2_voiced = mic2_aud(voiced_index+shift:voiced_index+shift+ms30-1);

[estimated_delay, GCCn,peak_position] = gcc_phat_alternative(mic1_voiced, mic2_voiced);
disp("Voiced:"+ estimated_delay/fs*1000);


% figure(1);
% subplot(311);
% plot(mic1_voiced);
% subtitle('Mic-1');
% subplot(312);
% plot(mic2_voiced);
% subtitle('Mic-2');
% subplot(3,1,3);
% plot(GCCn);
% hold on;
% scatter(peak_position, GCCn(peak_position), "filled");
% subtitle('Cross correlation function')


unvoiced_index = 22336;
mic1_unvoiced = mic1_aud(unvoiced_index:unvoiced_index+ms30-1);
mic2_unvoiced = mic2_aud(unvoiced_index+shift:unvoiced_index+shift+ms30-1);

[estimated_delay, GCCn, peak_position] = gcc_phat_alternative(mic1_unvoiced, mic2_unvoiced);
disp("Unoiced:"+ estimated_delay/fs*1000);

% 
% figure('DefaultAxesFontSize', 12);
% subplot(311);
% plot(mic1_unvoiced,'LineWidth',1.5);
% xlabel('Time (s)','FontSize', 12,'FontWeight','bold')
% ylabel('Amplitude','FontSize', 12,'FontWeight','bold')
% title('(a)','FontSize', 12,'FontWeight','bold');
% subplot(312);
% plot(mic2_unvoiced,'LineWidth',1.5);
% xlabel('Time (s)','FontSize', 12,'FontWeight','bold')
% ylabel('Amplitude','FontSize', 12,'FontWeight','bold')
% title('(b)','FontSize', 12,'FontWeight','bold');
% subplot(3,1,3);
% plot(GCCn,'LineWidth',1.5);
% hold on;
% scatter(peak_position, GCCn(peak_position), "filled");
% xlabel('Number of Samples','FontSize', 12,'FontWeight','bold')
% % ylabel('Amplitude','FontSize', 12,'FontWeight','bold')
% title('(c)','FontSize', 12,'FontWeight','bold');

n=300;
shift = 16;
ms30 = 0.03*fs;
audiolen = length(mic2_aud);
frameshift = 0.02*fs;

SNRs = [50,20,10,5,0,-5];
ste_array = csvread('csv_forPlots\ste_array_noNoise_2.csv');
ste_array=ste_array(:);

delay_array = zeros(length(SNRs), 160);

d_v = zeros(length(SNRs), 160);
d_uv = zeros(length(SNRs), 160);
for p=1:length(SNRs)

    SNR = SNRs(p);
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
        mic2_seg = mic2_aud(start2:end2);

        estimated_delay = gcc_phat_alternative(mic1_seg, mic2_seg);
        computed_delays(i) = estimated_delay;

        ste_frame = ste_array(i);
        if(ste_frame > 0.15)
            v(i) = estimated_delay;
        else
            u(i)=estimated_delay;
        end

    end
    delay_array(p,:) = computed_delays(1:160);
    d_v(p,:)=v(1:160);
    d_uv(p,:)=u(1:160);
end

% correct_estimates = sum(computed_delays(:) == shift);
% disp(correct_estimates);

frames_index = 1:160;

figure(1);
subplot(6,1,1);
scatter(frames_index,delay_array(1,:), 'filled');
ylim([-50 50])
subplot(612);
scatter(frames_index, delay_array(2,:), 'filled');
ylim([-50 50])
subplot(613);
scatter(frames_index, delay_array(3,:), 'filled');
ylim([-50 50])
subplot(614);
scatter(frames_index, delay_array(4,:), 'filled');
ylim([-50 50])
subplot(615);
scatter(frames_index, delay_array(5,:), 'filled');
ylim([-50 50])
subplot(616);
scatter(frames_index, delay_array(6,:), 'filled');
ylim([-50 50])


writematrix(delay_array, './csv_forPlots/phat_delay_array_arctic_b0399.csv');
writematrix(d_v, './csv_forPlots/phat_delay_vowel_arctic_b0399.csv');
% writematrix(delay_array, './csv_forPlots/phat_delay_array_2.csv');

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
% 


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
