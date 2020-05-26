clc; clear;
Eb_N0 = 0:40;
% ва 15/31
figure(1)
% ther_ber = berawgn(Eb_N0,'qam',4);
load('DataBase/CR=9_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'-',1.5,0,[0.45 0.45 0.45]);
load('DataBase/CR=15_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'-',1.5,0,[0.45 0.45 0.45]);
load('DataBase/CR=23_31_corM=1_2x2_RAYL_SPECIAL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=1');
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'-',1.5,0,[0.45 0.45 0.45]);
xlim([0 20 ])
legend ('MIMO 15/31','SISO15/31')
% legend ('MIMO 15/31','SISO15/31')
% figure(2)
% load('DataBase/CR=15_31_corM=1_2x2_RAYL_Wm=1_Ws=1_Mm=16_Ms=16_Exp=5.mat');
% plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0);
% plot_ber(ber_siso_mean,SNR,prm.bps_siso,'-',1.5,0,[0.45 0.45 0.45]);