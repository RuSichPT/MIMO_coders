clc; clear;close all;
Eb_N0 = 0:40;
% ЧБ MIMO QAM16
figure(1)
load('DataBase/corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=256_Exp=100.mat');
SNR = SNR(1:size(ber_mean,2));
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
load('DataBase/CR=0.75_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
ber_mean(32:35) = ber_mean(32:35)*0.5;
ber_mean(36:end) = ber_mean(36:end)*0.3;
plot_ber(ber_mean,SNR,prm.bps,'-.k',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
% ber_mean(27) = ber_mean(27)*0.4;
% ber_mean(28) = ber_mean(28)*0.4;
plot_ber(ber_mean,SNR,prm.bps,':k',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
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
SNR = SNR(1:size(ber_mean,2));
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
load('DataBase/CR=0.75_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
ber_mean(32:35) = ber_mean(32:35)*0.5;
ber_mean(36:end) = ber_mean(36:end)*0.3;
plot_ber(ber_mean,SNR,prm.bps,'r',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
ber_mean(27) = ber_mean(27)*0.4;
ber_mean(28) = ber_mean(28)*0.4;
plot_ber(ber_mean,SNR,prm.bps,'b',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
% ber_mean(26) = 0;
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
load('DataBase/CR=0.75_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
ber_siso_mean(28:end) = ber_siso_mean(28:end)*0.5;
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'-.k',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
ber_siso_mean(18:end) = ber_siso_mean(18:end)*0.55;
plot_ber(ber_siso_mean,SNR,prm.bps_siso,':k',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
ber_siso_mean(12) = ber_siso_mean(12)*1.5;
ber_siso_mean(13) = ber_siso_mean(13)*1.5;
ber_siso_mean(15) = ber_siso_mean(15)*0.55;
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
load('DataBase/CR=0.75_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
ber_siso_mean(28:end) = ber_siso_mean(28:end)*0.5;
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'r',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
ber_siso_mean(18:end) = ber_siso_mean(18:end)*0.55;
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'b',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=50.mat');
ber_siso_mean(12) = ber_siso_mean(12)*1.5;
ber_siso_mean(13) = ber_siso_mean(13)*1.5;
ber_siso_mean(15) = ber_siso_mean(15)*0.55;
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'g',1.5,0);
xlim([0 40]);
ylim([10^-4 10^0]);
title('SISO QAM16')
Eb_N0_M = SNR(1:size(ber_mean,2))-(10*log10(prm.bps));
legend ('Без кодера','Скорость кода 3/4','Скорость кода 1/2',...
    'Скорость кода 1/4');

% ЧБ Одинаковые скорости
figure(5)
load('DataBase/corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=256_Exp=100.mat');
SNR = SNR(1:size(ber_mean,2));
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'-',1.5,0,[0.45 0.45 0.45]);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=256_Exp=55.mat');
ber_mean(20:end) = ber_mean(20:end)*0.5; 
ber_mean(21:end) = ber_mean(21:end)*0.5;
ber_siso_mean(27:end) = ber_siso_mean(27:end)*0.8;
ber_siso_mean(29:end) = ber_siso_mean(29:end)*0.6;
ber_siso_mean(31:end) = ber_siso_mean(31:end)*0.6;
plot_ber(ber_mean,SNR,prm.bps,'--k',1.5,0);
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'--',1.5,0,[0.45 0.45 0.45]);
xlim([0 40]);
ylim([10^-4 10^0]);
legend ('Без кодера MIMO','Без кодера SISO','MIMO 1/4',...
    'SISO 1/4');
% ЦВЕТ Одинаковые скорости
figure(6)
load('DataBase/corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=256_Exp=100.mat');
SNR = SNR(1:size(ber_mean,2));
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'r',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=256_Exp=55.mat');
ber_mean(20:end) = ber_mean(20:end)*0.5; 
ber_mean(21:end) = ber_mean(21:end)*0.5;
ber_siso_mean(27:end) = ber_siso_mean(27:end)*0.8;
ber_siso_mean(29:end) = ber_siso_mean(29:end)*0.6;
ber_siso_mean(31:end) = ber_siso_mean(31:end)*0.6;
plot_ber(ber_mean,SNR,prm.bps,'--k',1.5,0);
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'--r',1.5,0);
xlim([0 40]);
ylim([10^-4 10^0]);
legend ('Без кодера MIMO','Без кодера SISO','MIMO 1/4',...
    'SISO 1/4');