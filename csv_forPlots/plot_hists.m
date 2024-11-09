clc;
clear all;
close all;


% gcc_d = csvread("csv_forPlots\phat_delay_array_arctic_b0399.csv");
% zff_d = csvread("csv_forPlots\zff_delay_array_arctic_b0399.csv");
% merged_d = csvread("csv_forPlots\merged_delay_array_arctic_b0399.csv");

gcc_d = csvread("csv_forPlots\phat_delay_vowel_arctic_b0399.csv");
zff_d = csvread("csv_forPlots\zff_delay_vowel_arctic_b0399.csv");
merged_d = csvread("csv_forPlots\merged_delay_vowel_arctic_b0399.csv");


nbins = 200;
gcc_20 = nonzeros(gcc_d(2,:));
zff_20 = nonzeros(zff_d(2,:));
merged_20 = nonzeros(merged_d(2,:));


% 20 dB hists
figure("Name","20dB");
% boxchart(delay_array');
subplot(311);
histogram(gcc_20,nbins);
subtitle('GCC');
xlim([-100 100]);
ylim([0 length(gcc_20)]);
subplot(312);
histogram(zff_20,nbins);
subtitle('ZFF');
xlim([-100 100]);
ylim([0 length(gcc_20)]);
subplot(313);
histogram(merged_20,nbins);
subtitle('Combined')
xlim([-100 100]);
ylim([0 length(gcc_20)]);


gcc_5 = nonzeros(gcc_d(4,:));
zff_5 = nonzeros(zff_d(4,:));
merged_5 = nonzeros(merged_d(4,:));

% 5 dB hists
figure(2);
% boxchart(delay_array');
subplot(311);
histogram(gcc_5,nbins);
subtitle('GCC');
xlim([-100 100]);
% ylim([0 length(gcc_20)]);
subplot(312);
histogram(zff_5,nbins);
subtitle('ZFF');
xlim([-100 100]);
% ylim([0 length(gcc_20)]);
subplot(313);
histogram(merged_5,nbins);
subtitle('Combined')
xlim([-100 100]);
% ylim([0 length(gcc_20)]);


gcc_neg5 = nonzeros(gcc_d(6,:));
zff_neg5 = nonzeros(zff_d(6,:));
merged_neg5 = nonzeros(merged_d(6,:));
% -5dB
figure(3);
% boxchart(delay_array');
subplot(311);
histogram(gcc_neg5,nbins);
subtitle('GCC');
xlim([-100 100]);
% ylim([0 length(gcc_20)]);
subplot(312);
histogram(zff_neg5,nbins);
subtitle('ZFF');
xlim([-100 100]);
% ylim([0 length(gcc_20)]);
subplot(313);
histogram(merged_neg5,nbins);
subtitle('Combined')
xlim([-100 100]);
% ylim([0 length(gcc_20)]);

