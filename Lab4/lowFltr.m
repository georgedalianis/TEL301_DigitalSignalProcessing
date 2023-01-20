function [h,fltr,w] = lowFltr(wc,fs,N,type)
%UNTITLED2 Summary of this function goes here
% types [rect hamm hann]
fc=wc/(2*pi);
Wn=fc/(fs/2);
if(type=='rect')
    fltr=fir1(N-1,Wn,rectwin(N));
elseif(type=='hamm')
    fltr=fir1(N-1,Wn,hamming(N));
elseif(type=='hann')
    fltr=fir1(N-1,Wn,hanning(N));
end
[h,w]=freqz(fltr,N);
end

