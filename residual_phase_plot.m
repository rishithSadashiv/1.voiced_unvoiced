%

clc;
clear all;
close all;




[d,fs]=audioread('./arctic_b0399.wav');

d=d(:,1);

noise = 0;
    
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

voiced_index = 5520;%11648;

mic1_voiced = mic1_aud(voiced_index:voiced_index+ms30-1);
mic2_voiced = mic2_aud(voiced_index+shift:voiced_index+shift+ms30-1);

[estimated_delay_v, GCCn_v,peak_position_v] = gcc_phat_alternative(mic1_voiced, mic2_voiced);
disp("Voiced:"+ estimated_delay_v/fs*1000);



unvoiced_index = 22336;
mic1_unvoiced = mic1_aud(unvoiced_index:unvoiced_index+ms30-1);
mic2_unvoiced = mic2_aud(unvoiced_index+shift:unvoiced_index+shift+ms30-1);

[estimated_delay_u, GCCn_u, peak_position_u] = gcc_phat_alternative(mic1_unvoiced, mic2_unvoiced);
disp("Unoiced:"+ estimated_delay_u/fs*1000);


res1=LPres(mic1_aud,fs,20,10,10,1);
hilb1=HilbertEnv(res1,fs,1);
res1=res1(:);
hilb1=hilb1(:);
residualPhase1 = res1./hilb1;

res2=LPres(mic2_aud,fs,20,10,10,1);
hilb2=HilbertEnv(res2,fs,1);
res2=res2(:);
hilb2=hilb2(:);
residualPhase2 = res2./hilb2;

mic1_rp_v = residualPhase1(voiced_index:voiced_index+ms30-1);
mic2_rp_v = residualPhase2(voiced_index+shift:voiced_index+ms30+shift-1);

[estimated_delay_r, GCCn_r, peak_position_r] = gcc_phat_alternative(mic1_rp_v, mic2_rp_v);
disp("residual phase:"+ estimated_delay_r/fs*1000);


t = (1:ms30)/fs;

f = figure('DefaultAxesFontSize', 16,'Position',[1 1 1500 800]);
subplot(231);
plot(t,mic1_unvoiced,'LineWidth',1.5);
xlabel('Time (s)','FontSize', 16,'FontWeight','bold')
ylabel('Amplitude','FontSize', 16,'FontWeight','bold')
title('(a)','FontSize', 20,'FontWeight','bold');
subplot(234);
plot(GCCn_u,'LineWidth',1.5);
hold on;
scatter(peak_position_u, GCCn_u(peak_position_u), "filled");
title('(d)','FontSize', 16,'FontWeight','bold');
xlabel('No. of Samples','FontSize', 16,'FontWeight','bold');
subplot(232);
plot(t,mic1_voiced,'LineWidth',1.5);
xlabel('Time (s)','FontSize', 16,'FontWeight','bold')
ylabel('Amplitude','FontSize', 16,'FontWeight','bold')
title('(b)','FontSize', 20,'FontWeight','bold');
subplot(235);
plot(GCCn_v,'LineWidth',1.5);
hold on;
scatter(peak_position_v, GCCn_v(peak_position_v), "filled");
title('(e)','FontSize', 20,'FontWeight','bold');
xlabel('No. of Samples','FontSize', 16,'FontWeight','bold');
subplot(2,3,3);
plot(t,mic1_rp_v,'LineWidth',1.5);
xlabel('Time (s)','FontSize', 16,'FontWeight','bold')
ylabel('Amplitude','FontSize', 16,'FontWeight','bold')
title('(c)','FontSize', 20,'FontWeight','bold');
subplot(236);
plot(GCCn_r,'LineWidth',1.5);
hold on;
scatter(peak_position_r, GCCn_r(peak_position_r), 'filled');
title('(f)','FontSize', 20,'FontWeight','bold');
xlabel('No. of Samples','FontSize', 16,'FontWeight','bold');

exportgraphics(f,'paper_plots/gcc-phat-comp.png','Resolution',300);


