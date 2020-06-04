clc; clear;close all;
Eb_N0 = 0:40;
% ЧБ MIMO QAM16
figure(1)
load('DataBase/corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=256_Exp=100.mat');
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
load('DataBase/CR=0.75_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=100.mat');
plot_ber(ber_mean,SNR,prm.bps,'-.k',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=43.mat');
ber_mean = mean(ber(1:42,:),1);
plot_ber(ber_mean,SNR,prm.bps,':k',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=100.mat');
ber_mean(26) = 0;
plot_ber(ber_mean,SNR,prm.bps,'--k',1.5,0);
xlim([0 40]);
ylim([10^-4 10^0]);
title('MIMO QAM16')
Eb_N0_M = SNR(1:size(ber_mean,2))-(10*log10(prm.bps));
legend ('Без кодера','Скорость кода 3/4','Скорость кода 1/2',...
    'Скорость кода 1/4');

%  ЦВЕТ MIMO QAM16
figure(2)
load('DataBase/corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=256_Exp=100.mat');
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
load('DataBase/CR=0.75_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=100.mat');
plot_ber(ber_mean,SNR,prm.bps,'r',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=43.mat');
ber_mean = mean(ber(1:42,:),1);
plot_ber(ber_mean,SNR,prm.bps,'b',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=100.mat');
ber_mean(26) = 0;
plot_ber(ber_mean,SNR,prm.bps,'g',1.5,0);
xlim([0 40]);
ylim([10^-4 10^0]);
title('MIMO QAM16')
Eb_N0_M = SNR(1:size(ber_mean,2))-(10*log10(prm.bps));
legend ('Без кодера','Скорость кода 3/4','Скорость кода 1/2',...
    'Скорость кода 1/4');

% ЧБ SISO QAM16
figure(3)
load('DataBase/corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=4_Ms=16_Exp=100.mat');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'k',1.5,0);
load('DataBase/CR=0.75_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=100.mat');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'-.k',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=43.mat');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,':k',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=100.mat');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'--k',1.5,0);
xlim([0 40]);
ylim([10^-4 10^0]);
title('SISO QAM16')
Eb_N0_M = SNR(1:size(ber_mean,2))-(10*log10(prm.bps));
legend ('Без кодера','Скорость кода 3/4','Скорость кода 1/2',...
    'Скорость кода 1/4');

% ЦВЕТ SISO QAM16
figure(4)
load('DataBase/corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=4_Ms=16_Exp=100.mat');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'k',1.5,0);
load('DataBase/CR=0.75_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=100.mat');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'r',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=43.mat');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'b',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=100.mat');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'g',1.5,0);
xlim([0 40]);
ylim([10^-4 10^0]);
title('SISO QAM16')
Eb_N0_M = SNR(1:size(ber_mean,2))-(10*log10(prm.bps));
legend ('Без кодера','Скорость кода 3/4','Скорость кода 1/2',...
    'Скорость кода 1/4');

figure(5)
load('DataBase/corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=256_Exp=100.mat');
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'k',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=256_Exp=100.mat');
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'k',1.5,0);
xlim([0 40]);
ylim([10^-4 10^0]);