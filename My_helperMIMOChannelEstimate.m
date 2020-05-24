function hD = My_helperMIMOChannelEstimate(rxData,ltfSC, prm)
% Estimate channel from the preamble signal data tones

[~, nltf, numRx] = size(rxData);
% nltf should be == numSTS

% Transmitted pilot mapping sequences
P = helperGetP(prm.numSTS);
  
Puse = P'; % Extract and conjugate the P matrix 

denom = nltf.*ltfSC;

hD = complex(zeros(prm.numSC,prm.numSTS,numRx));
for i = 1:numRx
    rxsym = squeeze(rxData(:,(1:nltf),i)); % Symbols per receive antenna
    for j = 1:prm.numSTS
        hD(:,j,i) = rxsym*Puse(:,j)./denom;
    end
end

end

% [EOF]