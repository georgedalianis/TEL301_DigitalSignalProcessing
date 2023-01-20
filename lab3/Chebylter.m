function [fl,W,numZ,denZ]  = Chebylter(N,delta_p,Wc,fs)
%% UNTITLED5 Summary of this function goes here
fc=Wc/(2*pi);% Wc/(2*pi);
Wk=fc/(fs/2);

samples=256;
[numZ,denZ] = cheby1(N,delta_p,Wk,'high');
figure
W=0:1/(samples-1):1; 
fl=freqz(numZ,denZ,samples);
title(['Chebyshev highpass with N=' num2str(N)])
end

