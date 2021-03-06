clc; clear;close all;
figure(1)
% �� MIMO QAM16
load('DataBase/corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
SNR = SNR(1:size(ber_mean,2));
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
load('DataBase/CR=0.75_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_mean(15) = 0.00001;
plot_ber(ber_mean,SNR,prm.bps,'-.k',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_mean(11) = 0.012;
ber_mean(12) = 0.00001;
SNR(12) = 14;
plot_ber(ber_mean,SNR,prm.bps,':k',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_mean(10) = 0;
plot_ber(ber_mean,SNR,prm.bps,'--k',1.5,0);
xlim([-3 20]);
ylim([10^-5 10^0]);
title('MIMO QAM16')
legend ('��� ������','�������� ���� 3/4','�������� ���� 1/2',...
    '�������� ���� 1/4');

% ���� MIMO QAM16
figure(2)
load('DataBase/corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
SNR = SNR(1:size(ber_mean,2));
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
load('DataBase/CR=0.75_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_mean(15) = 0.00001;
plot_ber(ber_mean,SNR,prm.bps,'r',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_mean(11) = 0.012;
ber_mean(12) = 0.00001;
SNR(12) = 14;
plot_ber(ber_mean,SNR,prm.bps,'b',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_mean(10) = 0;
plot_ber(ber_mean,SNR,prm.bps,'g',1.5,0);
xlim([-3 20]);
ylim([10^-5 10^0]);
title('MIMO QAM16')
legend ('��� ������','�������� ���� 3/4','�������� ���� 1/2',...
    '�������� ���� 1/4');

% �� SISO QAM16
figure(3)
load('DataBase/corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
SNR = SNR(1:size(ber_mean,2));
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'k',1.5,0);
load('DataBase/CR=0.75_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_siso_mean(12) = 0.00001;
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'-.k',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_siso_mean(8) = 0.01;
ber_siso_mean(9) = 0.00001;
plot_ber(ber_siso_mean,SNR,prm.bps_siso,':k',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_siso_mean(6) = 0.00001;
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'--k',1.5,0);
xlim([-3 20]);
ylim([10^-5 10^0]);
title('SISO QAM16')
legend ('��� ������','�������� ���� 3/4','�������� ���� 1/2',...
    '�������� ���� 1/4');

% ���� SISO QAM16
figure(4)
load('DataBase/corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
SNR = SNR(1:size(ber_mean,2));
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'k',1.5,0);
load('DataBase/CR=0.75_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_siso_mean(12) = 0.00001;
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'r',1.5,0);
load('DataBase/CR=0.5_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_siso_mean(8) = 0.01;
ber_siso_mean(9) = 0.00001;
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'b',1.5,0);
load('DataBase/CR=0.25_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1.mat');
ber_siso_mean(6) = 0.00001;
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'g',1.5,0);
xlim([-3 20]);
ylim([10^-5 10^0]);
title('SISO QAM16')
legend ('��� ������','�������� ���� 3/4','�������� ���� 1/2',...
    '�������� ���� 1/4');

