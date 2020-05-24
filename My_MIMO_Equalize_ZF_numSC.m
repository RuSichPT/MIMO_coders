function [Out_dataMod_ZF] = My_MIMO_Equalize_ZF_numSC(Y,H_estim)
% Модель Y = X*H+ksi; 

% Эквалайзер для каждой поднесущей
% Y - принятые символы [msc,symb_ofdm,numTx]
% H_estim - оценка матрицы вида [msc,numTx,numRx]
% msc - кол-во поднесущих,symb_ofdm - кол-во символов ofdm

for i = 1:size(Y,1)    
    h_estim = squeeze(H_estim(i,:,:));
    inv_H_ZF = inv(h_estim*h_estim');
    Out_dataMod_ZF(i,:,:) =  squeeze(Y(i,:,:))*h_estim'*inv_H_ZF;
end
% Выходные данные в символах IQ
end

