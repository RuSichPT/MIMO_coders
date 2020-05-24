function [H_WAV] = H_WAV_my_mimo(H_LS,varargin)
% varargin{1}-LEVEL
if nargin==2
    for i = 1:size(H_LS,2)
        H_WAV_RE = wdenoise(real(squeeze(H_LS(:,i,:))),varargin{1},'Wavelet','db4','DenoisingMethod','UniversalThreshold');
        H_WAV_IM = wdenoise(imag(squeeze(H_LS(:,i,:))),varargin{1},'Wavelet','db4','DenoisingMethod','UniversalThreshold');
        H_WAV(:,i,:) = H_WAV_RE+1i*H_WAV_IM;
    end
else
    N = length(H_LS);
    val_min = min(wmaxlev(N,'db4'),floor(log2(N)));
    for i = 1:size(H_LS,2)
        H_WAV_RE = wdenoise(real(squeeze(H_LS(:,i,:))),val_min,'Wavelet','db4','DenoisingMethod','UniversalThreshold');
        H_WAV_IM = wdenoise(imag(squeeze(H_LS(:,i,:))),val_min,'Wavelet','db4','DenoisingMethod','UniversalThreshold');
        H_WAV(:,i,:) = H_WAV_RE+1i*H_WAV_IM;
    end
end
end

