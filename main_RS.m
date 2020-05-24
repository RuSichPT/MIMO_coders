%% ---------Модель MIMO and SISO RS-------- 
clear;clc;close all;
%% Управление
flag_chanel = 'RAYL'; % 'AWGN' ,'RAYL','STATIC', 'BAD'
flag_cor_MIMO = 1; % 1-коррекция АЧХ (эквалайзер для MIMO) 2-Аламоути
flag_cor_SISO = 1; % коррекция АЧХ (эквалайзер для SISO)
flag_wav_MIMO = 0; % вейвлет шумоподавление для MIMO
flag_wav_SISO = 0; % вейвлет шумоподавление для SISO
flag_coder_RS = 1; % кодер RS вкл/выкл
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
M_RS = 5; % Кол-во бит на символ , задаем поле Галуа, степень полинома ДБ 3 <= M_RS <= 16.
N_RS_max = 2^M_RS-1;
K_RS = 16; % Кол-во символов RS сообщения ДБ N_RS-K_RS= 2t, где t - кол-во исправимых ошибок
N_RS = N_RS_max; % Кол-во символов RS кодового слова  ДБ N_RS_max>=N_RS 
%% Параметры OFDM 
prm.numSC = 234; % Кол-во поднессущих
prm.N_FFT = 256; % Длина FFT для OFDM
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
if flag_coder_RS == 1
    % Объекты кодера Рида-Соломона
    if ((N_RS-K_RS)<2) || (N_RS>N_RS_max)
        error('Проверь условия N_RS-K_RS > 2 and N_RS_max>N_RS');
    end
    rsEncoder = comm.RSEncoder('BitInput',true,'CodewordLength',N_RS,'MessageLength',K_RS);
    rsEncoder_siso = clone(rsEncoder);
    rsDecoder = comm.RSDecoder('BitInput',true,'CodewordLength',N_RS,'MessageLength',K_RS);
    rsDecoder_siso = clone(rsDecoder);
    prm.CodeRate = K_RS/N_RS;
    prm.Nsymb_ofdm = prm.Nsymb_ofdm/10*N_RS;
end
%------- 
prm.n = prm.bps*prm.Nsymb_ofdm*prm.numSC*numTx*prm.CodeRate;% Длина бинарного потока
prm.n_pilot = prm.Nsymb_ofdm_p*prm.numSC; % Кол-во бит на пилоты SISO
prm.n_siso = prm.bps_siso*prm.Nsymb_ofdm*prm.numSC*prm.CodeRate;% Длина бинарного потока
%% ---------Сам скрипт--------
SNR = 0:30;
Exp = 50;% Кол-во опытов
for indExp = 1:Exp
    %% Создание канала
    [H,H_siso] = create_chanel(flag_chanel,prm);  
    for indSNR = 1:length(SNR) 
        %% Формируем данные 
        if flag_coder_RS == 1
            Inp_data = randi([0 1],prm.n,1); % Передаваемые данные
            EncData = rsEncoder(Inp_data); % Кодер RS
        else
            Inp_data = randi([0 1],prm.n,1); % Передаваемые данные
            EncData = Inp_data;
        end
        Inp_Mat = reshape(EncData,length(EncData)/prm.bps,prm.bps); %Группируем биты 
        Inp_Sym = bi2de(Inp_Mat);  % Входные данные в символах
        % SISO
        if flag_coder_RS == 1
            Inp_data_siso = randi([0 1],prm.n_siso,1); % Передаваемые данные
            EncData_siso = rsEncoder_siso(Inp_data_siso); % Кодер RS
        else
            Inp_data_siso = randi([0 1],prm.n_siso,1); % Передаваемые данные
            EncData_siso = Inp_data_siso;
        end          
        Inp_Mat_siso = reshape(EncData_siso,length(EncData_siso)/prm.bps_siso,prm.bps_siso); %Группируем биты 
        Inp_Sym_siso = bi2de(Inp_Mat_siso);  % Входные данные в символах
        % Формируем пилоты для SISO 
        Inp_data_pilot = randi([0 1],prm.n_pilot,1);%  набор пилотов в битах
        %% Модулятор
        % MIMO
        Mod_data_inp_tmp = qammod(Inp_Sym,prm.M);% Модулятор QAM-M для полезной инф
        if flag_cor_MIMO == 2
            Mod_data_inp_tmp = ostbcEnc(Mod_data_inp_tmp);
        end 
        Mod_data_inp = reshape(Mod_data_inp_tmp,prm.numSC,prm.Nsymb_ofdm,prm.numTx);
        % Модулятор пилотов  MIMO
        [preambula,ltfSC] = My_helperGenPreamble(prm);
        % Модулятор пилотов  SISO
        Mod_data_inp_pilot = pskmod(Inp_data_pilot,2);
        % SISO
        Mod_data_inp_tmp_siso = qammod(Inp_Sym_siso,prm.M_siso);% Модулятор QAM-M для полезной инф
        Mod_data_inp_tmp_siso1 = cat(1,Mod_data_inp_pilot,Mod_data_inp_tmp_siso );
        Mod_data_inp_siso = reshape(Mod_data_inp_tmp_siso1,prm.numSC,prm.Nsymb_ofdm+prm.Nsymb_ofdm_p,1);
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
            Mod_data_out_ZF = reshape(Mod_data_out_ZF_tmp,prm.numSC*prm.Nsymb_ofdm,prm.numRx);
        elseif flag_cor_MIMO == 2
            Mod_data_out_ZF_tmp = Mod_data_out(:,prm.numTx+1:end,:);
            for i = 1:prm.Nsymb_ofdm
                Mod_data_out_ZF(:,i) = ostbcComb(squeeze(Mod_data_out_ZF_tmp(:,i,:)),H_estim);
            end
        else
            Mod_data_out_ZF_tmp = Mod_data_out(:,prm.numTx+1:end,:);
            Mod_data_out_ZF = reshape(Mod_data_out_ZF_tmp,prm.numSC*prm.Nsymb_ofdm,prm.numRx);
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
        Out_Sym = qamdemod(Mod_data_out_tmp,prm.M);    
        Mod_data_out_siso_tmp = Mod_data_out_siso_ZF(:);
        Out_Sym_siso = qamdemod(Mod_data_out_siso_tmp,prm.M_siso);
        %% Выходные данные
        Out_Mat = de2bi(Out_Sym);
        Out_data = Out_Mat(:);
        if(flag_coder_RS == 1)
            Out_data = rsDecoder(Out_data); % Декодер RS
        end
        Out_Mat_siso = de2bi(Out_Sym_siso);
        Out_data_siso  = Out_Mat_siso(:);
        if(flag_coder_RS == 1)
            Out_data_siso = rsDecoder_siso(Out_data_siso); % Декодер RS
        end
        
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
ther_ber_1 = berfading(Eb_N0,'qam',4,1);
ther_ber_4 = berfading(Eb_N0,'qam',4,4);
% ther_ber_1 = berawgn(Eb_N0,'qam',4);
figure(1)
plot_ber(ther_ber_1,Eb_N0,1,'g',1.5,0)
plot_ber(ther_ber_4,Eb_N0,1,'r',1.5,0)
plot_ber(ber_mean,SNR,prm.bps,'k',1.5,0)
plot_ber(ber_siso_mean,SNR,prm.bps_siso,'b',1.5,0)
% semilogy(ther_ber_4.data{1},ther_ber_4.data{2},'r','LineWidth',1.5)
legend("Теоретическая order = 1","Теоретическая order = 4","MIMO","SISO")%,"Теоретическая order = 4")