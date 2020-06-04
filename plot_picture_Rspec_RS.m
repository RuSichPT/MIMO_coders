clc; clear;%close all;
Eb_N0 = 0:40;
% ЧБ MIMO QAM16
figure(1)
load('DataBase/corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=30.mat');
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
load('DataBase/CR=23_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_mean,SNR,prm.bps,'-.k',1.5,0);
load('DataBase/CR=15_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_mean,SNR,prm.bps,':k',1.5,0);
load('DataBase/CR=9_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_mean,SNR,prm.bps,'--k',1.5,0);
ther_ber = berawgn(Eb_N0,'qam',16);
plot_ber(ther_ber,Eb_N0,1,'-',1.5,0,[0.45 0.45 0.45]);
xlim([0 22]);
ylim([10^-5 10^0]);
title('MIMO QAM16')
legend ('Без кодера','Скорость кода 3/4 (23/31)','Скорость кода 1/2 (15/31)',...
    'Скорость кода 1/4 (9/31)','Теоретическая awgn');

% ЦВЕТ MIMO QAM16
figure(2)
load('DataBase/corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=30.mat');
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
load('DataBase/CR=23_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_mean,SNR,prm.bps,'r',1.5,0);
load('DataBase/CR=15_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_mean,SNR,prm.bps,'b',1.5,0);
load('DataBase/CR=9_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_mean,SNR,prm.bps,'g',1.5,0);
ther_ber = berawgn(Eb_N0,'qam',16);
plot_ber(ther_ber,Eb_N0,1,'m',1.5,0);
xlim([0 22]);
ylim([10^-5 10^0]);
title('MIMO QAM16')
legend ('Без кодера','Скорость кода 3/4 (23/31)','Скорость кода 1/2 (15/31)',...
    'Скорость кода 1/4 (9/31)','Теоретическая awgn');

% ЧБ SISO QAM16
figure(3)
load('DataBase/corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=30.mat');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'k',1.5,0);
load('DataBase/CR=23_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'-.k',1.5,0);
load('DataBase/CR=15_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,':k',1.5,0);
load('DataBase/CR=9_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'--k',1.5,0);
ther_ber = berawgn(Eb_N0,'qam',16);
plot_ber(ther_ber,Eb_N0,1,'-',1.5,0,[0.45 0.45 0.45]);
xlim([0 22]);
ylim([10^-5 10^0]);
title('SISO QAM16')
legend ('Без кодера','Скорость кода 3/4 (23/31)','Скорость кода 1/2 (15/31)',...
    'Скорость кода 1/4 (9/31)','Теоретическая awgn');

% ЦВЕТ SISO QAM16
figure(4)
load('DataBase/corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=30.mat');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'k',1.5,0);
load('DataBase/CR=23_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'r',1.5,0);
load('DataBase/CR=15_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'b',1.5,0);
load('DataBase/CR=9_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'g',1.5,0);
ther_ber = berawgn(Eb_N0,'qam',16);
plot_ber(ther_ber,Eb_N0,1,'m',1.5,0);
xlim([0 22]);
ylim([10^-5 10^0]);
title('SISO QAM16')
legend ('Без кодера','Скорость кода 3/4 (23/31)','Скорость кода 1/2 (15/31)',...
    'Скорость кода 1/4 (9/31)','Теоретическая awgn');
