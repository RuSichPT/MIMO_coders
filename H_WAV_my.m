function [H_WAV] = H_WAV_my(H_LS,varargin)
% varargin{1}-LEVEL
if nargin==2
    H_WAV_RE = wdenoise(real(H_LS),varargin{1},'Wavelet','db4','DenoisingMethod','UniversalThreshold');
    H_WAV_IM = wdenoise(imag(H_LS),varargin{1},'Wavelet','db4','DenoisingMethod','UniversalThreshold');
    H_WAV=H_WAV_RE+1i*H_WAV_IM;
else
    N = length(H_LS);
    val_min = min(wmaxlev(N,'db4'),floor(log2(N)));
    H_WAV_RE = wdenoise(real(H_LS),val_min,'Wavelet','db4','DenoisingMethod','UniversalThreshold');
    H_WAV_IM = wdenoise(imag(H_LS),val_min,'Wavelet','db4','DenoisingMethod','UniversalThreshold');
    H_WAV=H_WAV_RE+1i*H_WAV_IM;
end

end

