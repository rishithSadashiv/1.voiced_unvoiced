%


clc;
clear all;
close all;


gcc_d = csvread("csv_forPlots\phat_delay_array.csv");
zff_d = csvread("csv_forPlots\zff_delay_array.csv");
merged_d = csvread("csv_forPlots\merged_delay_array.csv");

ste_array = csvread("csv_forPlots\ste_array_noNoise.csv");

SNRs = [50,20,10,5,0,-5];
frames_index = 1:160;
true_tdoa = 16.*ones(160,1);


f = figure('DefaultAxesFontSize', 20,'Position',[1 1 1500 800]);
subplot(331)    
scatter(frames_index,gcc_d(2,:), 12,'filled');
hold on;
plot(true_tdoa,"-.","Color","#00841a")
title("(a)","FontSize",20,"FontWeight","bold");
ylim([-50 50])
subplot(332);
scatter(frames_index,gcc_d(4,:), 12,'filled');
hold on;
plot(true_tdoa,"-.","Color","#00841a")
ylim([-50 50])
title("(b)","FontSize",20,"FontWeight","bold");
subplot(333);
scatter(frames_index,gcc_d(6,:), 12,'filled');
hold on;
plot(true_tdoa,"-.","Color","#00841a")
ylim([-50 50])
title("(c)","FontSize",20,"FontWeight","bold");
% subplot(4,3,10);
% plot(ste_array(1:160),'LineWidth',1.5);

subplot(334);
scatter(frames_index,zff_d(2,:), 12,'filled');
hold on;
plot(true_tdoa,"-.","Color","#00841a")
ylim([-50 50])
subplot(335);
scatter(frames_index,zff_d(4,:), 12,'filled');
hold on;
plot(true_tdoa,"-.","Color","#00841a")
ylim([-50 50])
subplot(336);
scatter(frames_index,zff_d(6,:), 12,'filled');
hold on;
plot(true_tdoa,"-.","Color","#00841a")
ylim([-50 50])
% subplot(4,3,11);
% plot(ste_array(1:160),'LineWidth',1.5);
% xlabel('Frame index')
subplot(337);
scatter(frames_index,merged_d(2,:), 12,'filled');
hold on;
plot(true_tdoa,"-.","Color","#00841a")
ylim([-50 50])
xlabel('Frame index','FontSize', 16,'FontWeight','bold')
subplot(338);
scatter(frames_index,merged_d(4,:), 12,'filled');
hold on;
plot(true_tdoa,"-.","Color","#00841a")
ylim([-50 50])
xlabel('Frame index','FontSize', 16,'FontWeight','bold')
subplot(339);
scatter(frames_index,merged_d(6,:), 12,'filled');
hold on;
plot(true_tdoa,"-.","Color","#00841a")
ylim([-50 50])
xlabel('Frame index','FontSize', 16,'FontWeight','bold')
% subplot(4,3,12);
% plot(ste_array(1:160),'LineWidth',1.5);
% xlabel('Frame index')
% 
% t = tiledlayout(1,1,'Padding','tight');

% t.OuterPosition = [0.25 0.25 3 3];


exportgraphics(f,'paper_plots/comp_fig.png','Resolution',300);



% 
% subplot(7,1,1);
% scatter(frames_index,delay_array(1,:), 16,'filled');
% ylim([-50 50])
% subplot(712);
% scatter(frames_index, delay_array(2,:), 16, 'filled');
% ylim([-50 50])
% subplot(713);
% scatter(frames_index, delay_array(3,:), 16, 'filled');
% ylim([-50 50])
% subplot(714);
% scatter(frames_index, delay_array(4,:), 16, 'filled');
% ylim([-50 50])
% subplot(715);
% scatter(frames_index, delay_array(5,:), 16,'filled');
% ylim([-50 50])
% subplot(716);
% scatter(frames_index, delay_array(6,:), 16,'filled');
% ylim([-50 50])
% subplot(717);
% plot(ste_array(1:160));




