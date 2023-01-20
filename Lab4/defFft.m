function [outfft,faxis] = defFft(fs,N,xt)
% nmin = ceil(min(t)/48); nmax = floor(max(t)/48); N = nmin:nmax;
outfft=fftshift(fft(xt,N)*(1/fs));
faxis=[-(fs/2):(fs/N):(fs/2)-(fs/N)];
plot(faxis,abs(outfft));hold on
%stem(faxis,abs(outfft));
title('Sampling with T_s='+fs)
grid
end

