clc; close all; clear all;
%% 1 Butterworth
%       data       %%%%%%%%%
fs = 10000; %Hz
passBand=3000; %Hz
delta_p = 3; %dB
stopband=4000;%Hz
% attenuation
delta_s=30; %dB

MyButtFilter(passBand,stopband,delta_p,delta_s,fs,1);
delta_s=50; %dB
MyButtFilter(passBand,stopband,delta_p,delta_s,fs,1);


%% 2 highpas

N= [2 16];
Wc=2;% rad/sec cutoff
Ts=0.2 %s 
fs=1/Ts;
delta_p=3; %dB

[fl(1,:),w(1,:),numZ,denZ]=Chebylter(N(1),delta_p,Wc,fs);
close
[fl(2,:),w(2,:),numZ,denZ]=Chebylter(N(2),delta_p,Wc,fs);
plot(w(1,:),mag2db(abs(fl(1,:))),'DisplayName','N=2');hold on
plot(w(2,:),mag2db(abs(fl(2,:))),'DisplayName','N=16')
title('2. Chebyshev highpass filters N=[2 16]')


%% 3
% ----3 a----%
% x(t)=1+cos(1000t) + cos(1600t) + cos(30000t)
f1=500/pi;
f2=8000/pi;
f3=15000/pi;
fs=10000;
Ts=1/fs;
N=500;
n=0:N-1;
x1= 1+ cos(2*pi*f1*n*Ts) + cos(2*pi*f2*n*Ts) + cos(2*pi*f3*n*Ts)
passBand=3000; %Hz
delta_p = 3; %dB
stopband=4000;%Hz
% attenuation
delta_s=30; %dB
[numZ,denZ]=MyButtFilter(passBand,stopband,delta_p,delta_s,fs,0);
figure;
subplot(3,1,1);
plot(n,x1,'displayname','signal');hold on 
title('The first 500 samples of x(t)=1+cos(1000t) + cos(1600t) + 3cos(30000t))');
fl=filter(numZ,denZ,x1);
plot(n,fl,'displayname','filtered','linewidth',2)
title('low pass filter');
subplot(3,1,2);
defFft(fs,N,x1)
title(' without low pass filter');

subplot(3,1,3);

defFft(fs,N,fl)
title(' with low pass filter');

delta_s=50; %dB
[numZ,denZ]=MyButtFilter(passBand,stopband,delta_p,delta_s,fs,0);
figure;
subplot(3,1,1);
plot(n,x1,'displayname','signal');hold on 
title('The first 500 samples of x(t)=1+cos(1000t) + cos(1600t) + 3cos(30000t))');
fl=filter(numZ,denZ,x1);
plot(n,fl,'displayname','filtered','linewidth',2)
title('low pass filter');
subplot(3,1,2);
defFft(fs,N,x1)
title(' without low pass filter');

subplot(3,1,3);

defFft(fs,N,fl)
title(' with low pass filter');

%-----3b-----%
% x(t)=1+cos(1.5t)+cos(5t)
f1=0.75/pi;
f2=2.5/pi;
fs=5;
Ts=1/fs;
N=500;
n=0:N-1;
x2= 1+ cos(2*pi*f1*n*Ts) + cos(2*pi*f2*n*Ts)
% figure;
% subplot(3,1,1);
% stem(n,x2);
% title('The first 500 samples of x(t)=1+cos(1.5t)+cos(5t)');
% Chebylter(N,delta_p,Wc,fs)

N= [2 16];
Wc=2;% rad/sec cutoff
Ts=0.2 %s 
fs=1/Ts;
delta_p=3; %dB
clear fl w
for k=1:2
    
    [fl(k,:),w(k,:),numZ,denZ]=Chebylter(N(k),delta_p,Wc,fs);
    close
    figure
    subplot(3,1,1)
    plot(n,x2,'DisplayName','Original Signal');hold on
    %plot(w(k,:),mag2db(abs(fl(k,:))),'DisplayName',['N=' num2str(N(k))]);
    fl2=filter(numZ,denZ,x2);
    plot(n,fl2,'displayname','filtered','linewidth',2);
    title([' Chebyshev highpass filter N=' num2str(N(k))])
    subplot(3,1,2)
    defFft(fs,500,x2)
    
    
    title('low pass filter');
    subplot(3,1,3);
    defFft(fs,500,fl2)
    title(' without low pass filter');

    
    
    
end
