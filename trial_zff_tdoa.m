%


clc;
clear all;
close all;

[d,fs]=audioread('./arctic_b0399.wav');
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


ms30 = 0.03 * fs;
shift = 16;
voiced_index = 5520;%11648;
unvoiced_index = 22336;
t = (1:ms30)/fs;


mic1_voiced = mic1_aud(voiced_index:voiced_index+ms30-1);
mic2_voiced = mic2_aud(voiced_index+shift:voiced_index+shift+ms30-1);
% [m1_v_epochlocs,m1_zsp1,m1_vgclocssp1] = EpochsbyZFF(mic1_voiced, fs);
% m1_v_imp_train = zeros(1,ms30);
% m1_v_imp_train(m1_v_epochlocs1) = 1;

% disp(length(m1_imp_train));
% 
% [m2_v_epochlocs,m2_zsp1,m2_vgclocssp1] = EpochsbyZFF(mic2_voiced, fs);
% m2_imp_train = zeros(1,ms30);
% m2_imp_train(m2_v_epochlocs) = 1;
[m1_voiced_spec,m1_f] = pspectrum(mic1_voiced, fs);
[m2_voiced_spec,m2_f] = pspectrum(mic2_voiced, fs);


unvoiced_index = 22336;
mic1_unvoiced = mic1_aud(unvoiced_index:unvoiced_index+ms30-1);
mic2_unvoiced = mic2_aud(unvoiced_index+shift:unvoiced_index+shift+ms30-1);
% [m1_un_epochlocs,m1_un_zsp1,m1_un_vgclocssp1] = EpochsbyZFF(mic1_unvoiced, fs);
% m1_un_imp_train = zeros(1,ms30);
% m1_un_imp_train(m1_un_epochlocs) = 1;
% 
% [m2_un_epochlocs,m2_un_zsp1,m2_un_vgclocssp1] = EpochsbyZFF(mic2_unvoiced, fs);
% m2_un_imp_train = zeros(1,ms30);
% m2_un_imp_train(m2_un_epochlocs) = 1;
[m1_un_voiced_spec,m1_f] = pspectrum(mic1_unvoiced, fs);
[m2_un_voiced_spec,m2_f] = pspectrum(mic2_unvoiced, fs);


[m1_epochlocs,m1_zsp1,m1_vgclocssp1] = EpochsbyZFF(mic1_aud, fs);
m1_imp_train = zeros(1,length(mic1_aud));
m1_imp_train(m1_epochlocs) = 1;
m1_imp_train_v = m1_imp_train(voiced_index:voiced_index+ms30-1);
m1_zf_v = m1_zsp1(voiced_index:voiced_index+ms30-1);
m1_imp_train_uv = m1_imp_train(unvoiced_index:unvoiced_index+ms30-1);

[m2_epochlocs,m2_zsp1,m2_vgclocssp1] = EpochsbyZFF(mic2_aud, fs);
m2_imp_train = zeros(1,length(mic2_aud));
m2_imp_train(m2_epochlocs) = 1;
m2_imp_train_v = m2_imp_train(voiced_index+shift:voiced_index+shift+ms30-1);
m2_zf_v = m2_zsp1(voiced_index+shift:voiced_index+shift+ms30-1); 
m2_imp_train_uv = m2_imp_train(unvoiced_index+shift:unvoiced_index+shift+ms30-1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots
% figure(1);
% subplot(211);
% plot(mic1_voiced);
% title('Voiced segment: Mic 1');
% subplot(612);
% plot(m1_imp_train);
% title('ZFF impulse train: Mic 1');
% subplot(613);
% plot(mic2_voiced);
% title('Voiced segment: Mic 2');
% subplot(614);
% plot(m2_imp_train);
% title('ZFF impulse train: Mic 2');
% subplot(615);
% plot(m1_f,pow2db(m2_voiced_spec));
% xlabel('Frequency (Hz)');
% ylabel('Power Spectrum (dB)');
% title('Log Power Spectrum: Mic 1');
% subplot(212);
% plot(m1_f,pow2db(m1_voiced_spec));
% xlabel('Frequency (Hz)');
% ylabel('Power Spectrum (dB)');
% title('Log Power Spectrum: Mic 1');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots
% 
% figure(2);
% subplot(421);
% plot(mic1_voiced);
% title('Voiced Segment');
% subtitle('Mic 1');
% subplot(422);
% plot(mic1_unvoiced);
% title('Unvoiced Segment');
% subtitle('Mic 1');
% subplot(423);
% plot(m1_f, pow2db(m1_voiced_spec));
% subtitle('Log power spectrum: Mic 1')
% xlabel('Frequency (Hz)');
% ylabel('Power Spectrum (dB)');
% subplot(424);
% plot(m1_f, pow2db(m1_un_voiced_spec));
% subtitle('Log power spectrum: Mic 1');
% xlabel('Frequency (Hz)');
% ylabel('Power Spectrum (dB)');
% subplot(425);
% plot(mic2_voiced);
% subtitle('Mic 2');
% subplot(426);
% plot(mic2_unvoiced);
% subtitle('Mic 2');
% subplot(427);
% plot(m2_f, pow2db(m2_voiced_spec));
% subtitle('Log power spectrum: Mic 2');
% xlabel('Frequency (Hz)');
% ylabel('Power Spectrum (dB)');
% subplot(428);
% plot(m2_f, pow2db(m2_un_voiced_spec));
% subtitle('Log power spectrum: Mic 2');
% xlabel('Frequency (Hz)');
% ylabel('Power Spectrum (dB)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 1 %% used
% f = figure('DefaultAxesFontSize', 16,'Position',[1 1 1500 800]); 
% subplot(221);
% plot(t, mic1_voiced,'LineWidth',1.5);
% title('(a)', 'FontSize', 20,'FontWeight','bold');
% xlabel('Time (s)', 'FontSize', 16,'FontWeight','bold');
% ylabel('Amplitude', 'FontSize', 16,'FontWeight','bold');
% subplot(223);
% plot(t,mic1_unvoiced,'LineWidth',1.5);
% xlabel('Time (s)', 'FontSize', 16,'FontWeight','bold');
% ylabel('Amplitude', 'FontSize', 16,'FontWeight','bold');
% % title('Unvoiced Segment');
% title('(c)', 'FontSize', 20,'FontWeight','bold');
% subplot(222);
% plot(m1_f, pow2db(m1_voiced_spec),'LineWidth',1.5);
% % subtitle('Log power spectrum: Mic 1')
% xlabel('Frequency (Hz)', 'FontSize', 16,'FontWeight','bold');
% ylabel('Power Spectrum (dB)', 'FontSize', 16,'FontWeight','bold');
% title('(b)', 'FontSize', 20,'FontWeight','bold')
% subplot(224);
% plot(m1_f, pow2db(m1_un_voiced_spec),'LineWidth',1.5);
% title('(d)', 'FontSize', 20,'FontWeight','bold');
% xlabel('Frequency (Hz)', 'FontSize', 16,'FontWeight','bold');
% ylabel('Power Spectrum (dB)', 'FontSize', 16,'FontWeight','bold');
% 
% exportgraphics(f,'paper_plots/powerspec.png','Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot 1 %% used 

% 
% 
% figure(3);
% subplot(421);
% plot(mic1_voiced);
% title('Voiced Segment');
% subtitle('Mic 1');
% subplot(422);
% plot(mic1_unvoiced);
% title('Unvoiced Segment');
% subtitle('Mic 1');
% subplot(423);
% plot(m1_imp_train);
% subtitle('ZFF impulse train: Mic 1')
% subplot(424);
% plot(m1_un_imp_train);
% subtitle('ZFF impulse train: Mic 1');
% subplot(425);
% plot(mic2_voiced);
% subtitle('Mic 2');
% subplot(426);
% plot(mic2_unvoiced);
% subtitle('Mic 2');
% subplot(427);
% plot(m2_imp_train);
% subtitle('ZFF impulse train: Mic 2');
% subplot(428);
% plot(m2_un_imp_train);
% subtitle('ZFF impulse train: Mic 2');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot 2 %% used
i1 = find(m1_imp_train_v);
i2 = find(m2_imp_train_v);

f=figure('DefaultAxesFontSize', 12,'Position',[1 1 1500 800]);
subplot(421);
plot(t,mic1_voiced,'LineWidth',1.5);
xlabel('Time (s)','FontSize', 12,'FontWeight','bold')
ylabel('Amplitude','FontSize', 12,'FontWeight','bold')
title('(a)','FontSize', 12,'FontWeight','bold');
subplot(422);
plot(t,mic2_voiced,'LineWidth',1.5);
xlabel('Time (s)','FontSize', 12,'FontWeight','bold')
ylabel('Amplitude','FontSize', 12,'FontWeight','bold')
title('(b)','FontSize', 12,'FontWeight','bold');
subplot(4,2,3);
plot(t,m1_zf_v,'LineWidth',1.5);
hold on;
xline(i1/fs, '--b','LineWidth',1.5);
xlabel('Time (s)','FontSize', 12,'FontWeight','bold')
ylabel('Amplitude','FontSize', 12,'FontWeight','bold')
title('(c)','FontSize', 12,'FontWeight','bold');
subplot(424);
plot(t,m2_zf_v,'LineWidth',1.5);
hold on;
xline(i2/fs, '--r','LineWidth',1.5);
xlabel('Time (s)','FontSize', 12,'FontWeight','bold')
ylabel('Amplitude','FontSize', 12,'FontWeight','bold')
title('(d)','FontSize', 12,'FontWeight','bold');
subplot(425);
plot(t,m1_imp_train_v,'LineWidth',1.5);
xlabel('Time (s)','FontSize', 12,'FontWeight','bold')
ylabel('Amplitude','FontSize', 12,'FontWeight','bold')
title('(e)','FontSize', 12,'FontWeight','bold');
subplot(426);
plot(t,m2_imp_train_v, 'red','LineWidth',1.5);
% xline(i2/fs, '--r','LineWidth',1.5);
xlabel('Time (s)','FontSize', 12,'FontWeight','bold')
ylabel('Amplitude','FontSize', 12,'FontWeight','bold')
title('(f)','FontSize', 12,'FontWeight','bold');
subplot(4,2,7:8);
plot(t,m1_imp_train_v,'LineWidth',1.5);
hold on;
plot(t,m2_imp_train_v, 'red','LineWidth',1.5);
% xline(i2/fs, '--r','LineWidth',1.5);
xlabel('Time (s)','FontSize', 12,'FontWeight','bold')
ylabel('Amplitude','FontSize', 12,'FontWeight','bold')
title('(g)','FontSize', 12,'FontWeight','bold');

exportgraphics(f,'paper_plots/zff_plots.png','Resolution',300);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot 2 %% used


% figure(6);
% subplot(411);
% plot(mic1_voiced);
% subplot(412);
% plot(mic2_voiced);
% subplot(413);
% plot(mic1_unvoiced);
% subplot(414);
% plot(mic2_unvoiced);

% [a, b] = xcorr(m1_imp_train_v, m2_imp_train_v);
% figure(5);
% plot(b,a);

% [value, peak_positions] = max(a);
% peak_position = max(peak_positions);
% delay = b(peak_position);
% disp("With xcorr function:"+delay);


delay_v = computeDelayManual(m1_imp_train_v,m2_imp_train_v);
disp("With manual method (voiced):"+delay_v/fs*1000);
delay_uv = computeDelayManual(m1_imp_train_uv,m2_imp_train_uv);
disp("With manual method (unvoiced):"+delay_uv/fs*1000);



% https://www.scicoding.com/convolution-in-python-3-essential-packages/
% https://in.mathworks.com/matlabcentral/answers/365308-manual-code-for-convolution
% x=m1_imp_train_v;
% h=m2_imp_train_v;
% m=length(x);
% n=length(h);
% l = n+m-1;
% X=[x,zeros(1,l-n)]; 
% H=[h,zeros(1,l-m)]; 
% for i=1:n+m-1
% Y(i)=0;
% for j=1:m
% if(i-j+1>0)
% Y(i)=Y(i)+abs(X(j)-H(i-j+1));
% else
% end
% end
% end
% Y=Y;

% figure(10);
% plot(Y)



% function [estimated_delay] = meanAbsDev(imp1, imp2)
% 
% end

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

