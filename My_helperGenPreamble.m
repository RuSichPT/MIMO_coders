function [y,ltfSC] = My_helperGenPreamble(prm)
% Generate the Preamble signal for channel estimation.
Nltf = prm.numSTS; % number of preamble symbols

% Frequency subcarrier tones
x = randi([0 1],prm.numSC,1);
ltfSC = qammod(x,2);

P = helperGetP(prm.numSTS);    
Pred = P;

% Define LTF(Long training field) and output variable sizes
ltfTx = complex(zeros(prm.N_FFT,prm.numSTS));
symLen = prm.N_FFT+prm.CyclicPrefixLength;

% Generate and modulate each LTF symbol
y = complex(zeros(symLen*Nltf,prm.numSTS));
for i = 1:Nltf  
    ltfTx = ltfSC*Pred(:, i).';%(prm.CarriersLocations,:)
    % OFDM modulation
    tmp = ofdmmod(reshape(ltfTx, [prm.numSC,1,prm.numSTS]), ...
        prm.N_FFT, prm.CyclicPrefixLength,prm.NullCarrierIndices);
    y((i-1)*symLen+(1:symLen),:) = tmp;
end
