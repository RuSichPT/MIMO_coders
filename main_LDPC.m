%% ---------Модель MIMO and SISO LDPC-------- 
clear;clc;close all;
%% Управление
flag_chanel = 'RAYL'; % 'AWGN' ,'RAYL','STATIC', 'BAD'
flag_cor_MIMO = 2; % 1-коррекция АЧХ (эквалайзер для MIMO) 2-Аламоути
flag_cor_SISO = 1; % коррекция АЧХ (эквалайзер для SISO)
flag_wav_MIMO = 0; % вейвлет шумоподавление для MIMO
flag_wav_SISO = 0; % вейвлет шумоподавление для SISO
flag_coder_LDPC = 1; % кодер LDPC вкл/выкл
%% Параметры системы MIMO
prm.numTx = 2; % Кол-во излучающих антен
prm.numRx = 2; % Кол-во приемных антен
prm.numSTS = prm.numTx; % Кол-во потоков
prm.M = 4;% Порядок модуляции
prm.bps = log2(prm.M); % Коль-во бит на символ в секунду
prm.CodeRate = 1; % CONST не менять. по умолчанию
%% Параметры системы SISO
prm.M_siso = 4;% Порядок модуляции
prm.bps_siso = log2(prm.M_siso); % Коль-во бит на символ в секунду
prm.Nsymb_ofdm_p = 1; % Кол-во пилотных символов OFDM 
%% Параметры кодера
codeRate = 1/2; % 1/2, 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 5/6, 8/9, or 9/10
CODEWORD_LENGTH = 64800;
N_CodeWord = 2; % Кол-во кодовых слов
ofdm_symb = 72; % ДБ n = 64800; 
%% Параметры OFDM 
prm.numSC = 450; % Кол-во поднессущих
prm.N_FFT = 512; % Длина FFT для OFDM
prm.Nsymb_ofdm = 100; % Кол-во символов OFDM от каждой антенны
prm.CyclicPrefixLength = 64;  % длина защитных интервалов = 2*Ngi
prm.tmp_NCI = prm.N_FFT - prm.numSC;
prm.NullCarrierIndices = [1:prm.tmp_NCI/2 prm.N_FFT-prm.tmp_NCI/2+1:prm.N_FFT]'; % Guards and DC
%% Параметры канала
prm.SampleRate = 40e6;
dt = 1/prm.SampleRate;
prm.tau_mimo = [0 1*dt 2*dt];
prm.pdB_mimo = [0 -9 -12];
prm.tau_siso = [0 1*dt 2*dt];
prm.pdB_siso = [0 -9 -12];
% prm.tau_mimo = 0;
% prm.pdB_mimo = 0;
% prm.tau_siso = 0;
% prm.pdB_siso = 0;
if flag_cor_MIMO == 2
    ostbcEnc = comm.OSTBCEncoder('NumTransmitAntennas',prm.numTx);
    ostbcComb = comm.OSTBCCombiner('NumReceiveAntennas',prm.numRx);
    numTx = 1;
else
    numTx = prm.numTx;
end
if flag_coder_LDPC == 1
    % Объекты кодера LDPC
    LDPCParityCheckMatrix = dvbs2ldpc(codeRate);
    ldpcEnc = comm.LDPCEncoder();
    ldpcDec = comm.LDPCDecoder();
    prm.Nsymb_ofdm = ofdm_symb;
    prm.CodeRate = codeRate;
    Nsymb_ofdm_mimo = ofdm_symb/numTx;
    prm.Nsymb_ofdm = ofdm_symb;
else
    Nsymb_ofdm_mimo = prm.Nsymb_ofdm;
    Nsymb_ofdm_siso = prm.Nsymb_ofdm;
    N_CodeWord = 1;
end
%------- 
prm.n = prm.bps*Nsymb_ofdm_mimo*prm.numSC*numTx*prm.CodeRate;% Длина бинарного потока
prm.n_pilot = prm.Nsymb_ofdm_p*prm.numSC; % Кол-во бит на пилоты SISO
prm.n_siso = prm.bps_siso*prm.Nsymb_ofdm*prm.numSC*prm.CodeRate;% Длина бинарного потока
%% ---------Сам скрипт--------
SNR = 0:30;
Exp = 1;% Кол-во опытов
for indExp = 1:Exp
    %% Создание канала
    [H,H_siso] = create_chanel(flag_chanel,prm);  
    for indSNR = 1:length(SNR) 
        %% Формируем данные 
        if flag_coder_LDPC == 1
            Inp_data = randi([0 1],prm.n*N_CodeWord,1); % Передаваемые данные
            EncData = [];
            for i=1:N_CodeWord
                EncData = cat(1,EncData,ldpcEnc(Inp_data(1+(i-1)*prm.n:prm.n+(i-1)*prm.n))); % Кодер LDPC
            end
        else
            Inp_data = randi([0 1],prm.n,1); % Передаваемые данные
            EncData = Inp_data;
        end
        % SISO
        if flag_coder_LDPC == 1
            Inp_data_siso = randi([0 1],prm.n_siso*N_CodeWord,1); % Передаваемые данные
            EncData_siso = [];
            for i=1:N_CodeWord
                EncData_siso = cat(1,EncData_siso,ldpcEnc(Inp_data_siso(1+(i-1)*prm.n:prm.n+(i-1)*prm.n))); % Кодер LDPC
            end           
        else
            Inp_data_siso = randi([0 1],prm.n_siso,1); % Передаваемые данные
            EncData_siso = Inp_data_siso;
        end          
        % Формируем пилоты для SISO 
        Inp_data_pilot = randi([0 1],prm.n_pilot,1);%  набор пилотов в битах
        %% Модулятор
        % MIMO
        Mod_data_inp_tmp = qammod(EncData,prm.M,'InputType','bit');% Модулятор QAM-M для полезной инф
        if flag_cor_MIMO == 2
            Mod_data_inp_tmp = ostbcEnc(Mod_data_inp_tmp);
        end 
        Mod_data_inp = reshape(Mod_data_inp_tmp,prm.numSC,Nsymb_ofdm_mimo*N_CodeWord,prm.numTx);
        % Модулятор пилотов  MIMO
        [preambula,ltfSC] = My_helperGenPreamble(prm);
        % Модулятор пилотов  SISO
        Mod_data_inp_pilot = pskmod(Inp_data_pilot,2);
        % SISO
        Mod_data_inp_tmp_siso = qammod(EncData_siso,prm.M_siso,'InputType','bit');% Модулятор QAM-M для полезной инф
        Mod_data_inp_tmp_siso1 = cat(1,Mod_data_inp_pilot,Mod_data_inp_tmp_siso );
        Mod_data_inp_siso = reshape(Mod_data_inp_tmp_siso1,prm.numSC,prm.Nsymb_ofdm*N_CodeWord+prm.Nsymb_ofdm_p,1);
        %% Модулятор OFDM
        OFDM_data = ofdmmod(Mod_data_inp,prm.N_FFT,prm.CyclicPrefixLength,...
                     prm.NullCarrierIndices);               
        OFDM_data_siso = ofdmmod(Mod_data_inp_siso,prm.N_FFT,prm.CyclicPrefixLength,...
                     prm.NullCarrierIndices);        
        OFDM_data = [preambula ; OFDM_data];        
        %% Прохождение канала
        switch flag_chanel
            case 'RAYL'
%                 H.Visualization = 'Impulse and frequency responses';
%                 H.AntennaPairsToDisplay = [2,1];
                [Chanel_data, H_ist] = H(OFDM_data);
                [Chanel_data_siso, H_ist_siso] = H_siso(OFDM_data_siso);            
            otherwise                  
                Chanel_data  = OFDM_data*H;
                Chanel_data_siso  = OFDM_data_siso*H_siso; 

        end   
        %% Собственный  шум
        [Noise_data,sigma] = my_awgn(Chanel_data,SNR(indSNR));%SNR(indSNR)
        [Noise_data_siso,sigma_siso] = my_awgn(Chanel_data_siso,SNR(indSNR));%SNR(indSNR)
        %% Демодулятор OFDM
        Mod_data_out = ofdmdemod(Noise_data,prm.N_FFT,prm.CyclicPrefixLength,prm.CyclicPrefixLength, ...
            prm.NullCarrierIndices);
        Mod_data_out_siso = ofdmdemod(Noise_data_siso,prm.N_FFT,prm.CyclicPrefixLength,prm.CyclicPrefixLength, ...
            prm.NullCarrierIndices);  
        %% Оценка канала  
        H_estim = My_helperMIMOChannelEstimate(Mod_data_out(:,1:prm.numSTS,:),ltfSC,prm);
        H_estim_siso = Mod_data_out_siso(:,1)./Mod_data_inp_pilot;
        Mod_data_out_siso(:,1:prm.Nsymb_ofdm_p) = [];
        %% Вейвлет шумоподавление
        if flag_wav_SISO == 1
            H_estim_siso = H_WAV_my(H_estim_siso);
        end
        if flag_wav_MIMO == 1
            H_estim = H_WAV_my_mimo(H_estim);
        end
        %% Эквалайзер 
        %ZF MIMO
        if flag_cor_MIMO == 1
            Mod_data_out_ZF_tmp= My_MIMO_Equalize_ZF_numSC(Mod_data_out(:,prm.numTx+1:end,:),H_estim);
            Mod_data_out_ZF = reshape(Mod_data_out_ZF_tmp,prm.numSC*Nsymb_ofdm_mimo*N_CodeWord,prm.numRx);
        elseif flag_cor_MIMO == 2
            Mod_data_out_ZF_tmp = Mod_data_out(:,prm.numTx+1:end,:);
            for i = 1:prm.Nsymb_ofdm*N_CodeWord
                Mod_data_out_ZF(:,i) = ostbcComb(squeeze(Mod_data_out_ZF_tmp(:,i,:)),H_estim);
            end
        else
            Mod_data_out_ZF_tmp = Mod_data_out(:,prm.numTx+1:end,:);
            Mod_data_out_ZF = reshape(Mod_data_out_ZF_tmp,prm.numSC*prm.Nsymb_ofdm*N_CodeWord,prm.numRx);
        end
        %ZF SISO
        if flag_cor_SISO == 1
            for i = 1:size(Mod_data_out_siso,2)
                Mod_data_out_siso_ZF(:,i) = Mod_data_out_siso(:,i)./H_estim_siso;
            end
        else
            Mod_data_out_siso_ZF = Mod_data_out_siso;
        end        
        %% Демодулятор
        Mod_data_out_tmp = Mod_data_out_ZF(:);  
        Mod_data_out_siso_tmp = Mod_data_out_siso_ZF(:);     
        if(flag_coder_LDPC == 1)
            CWL = CODEWORD_LENGTH;
            Out_encode = qamdemod(Mod_data_out_tmp,prm.M,'OutputType','approxllr');
            Out_data = [];
            for i = 1:N_CodeWord
                Out_data = cat(1,Out_data,ldpcDec(Out_encode(1+(i-1)*CWL:CWL+(i-1)*CWL)));% Декодер LDPC
            end
        else
            Out_data = qamdemod(Mod_data_out_tmp,prm.M,'OutputType','bit');  
        end
        if(flag_coder_LDPC == 1)
            Out_encode_siso = qamdemod(Mod_data_out_siso_tmp,prm.M_siso,'OutputType','approxllr');
            Out_data_siso = [];
            for i = 1:N_CodeWord
                Out_data_siso = cat(1,Out_data_siso,ldpcDec(Out_encode_siso(1+(i-1)*CWL:CWL+(i-1)*CWL))); % Декодер LDPC
            end
        else
            Out_data_siso = qamdemod(Mod_data_out_siso_tmp,prm.M_siso,'OutputType','bit');  
        end
        %% BER
        [~,ber(indExp,indSNR)] = biterr(Out_data,Inp_data);
        [~,ber_siso(indExp,indSNR)] = biterr(Out_data_siso,Inp_data_siso);
        fprintf('Complete %d db \n',SNR(indSNR)); 
    end
    fprintf('Complete %d  \n',indExp);
end
if flag_cor_MIMO ~= 2 
    prm.bps = prm.bps*prm.numTx;
end
ber_mean = mean(ber,1);
ber_siso_mean = mean(ber_siso,1);
Eb_N0 = 0:40;
% ther_ber_1 = berfading(Eb_N0,'qam',4,1);
ther_ber_4 = berfading(Eb_N0,'qam',4,4);
ther_ber_1 = berawgn(Eb_N0,'qam',4);
figure(1)
plot_ber(ther_ber_1,Eb_N0,1,'g',1.5,0)
plot_ber(ther_ber_4,Eb_N0,1,'r',1.5,0)
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0)
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'b',1.5,0)
% semilogy(ther_ber_4.data{1},ther_ber_4.data{2},'r','LineWidth',1.5)
legend("Теоретическая order = 1","Теоретическая order = 4","MIMO","SISO")%,"Теоретическая order = 4")