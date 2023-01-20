clc;clear all; close all;

%%  1. Lowpass rectangular and hamming 
wc=0.48*pi; fs = 100; N=21;
[h(1,:),fltr(1,:),w(1,:)]=lowFltr(wc,fs,N,'rect');
[h(2,:),fltr(2,:),w(2,:)]=lowFltr(wc,fs,N,'hamm');

figure
plot(w(1,:),h(1,:),'displayname','Rectagular');hold on
plot(w(2,:),h(2,:),'displayname','Hamming')
title('Rectagular and Hamming')
legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  2. Lowpass Hamming and Hanning 
wc=0.5*pi; fs = 100; N=[21 41];

%%% Hamming and Hanning     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for k=1:length(N)
    subplot(1,2,1)
        [h,f,w]=lowFltr(wc,fs,N(k),'hamm');
        plot(w,h,'displayname',['N=' int2str(N(k))]);hold on
        legend
        title('Hamming')
    subplot(1,2,2)
        [h,f,w]=lowFltr(wc,fs,N(k),'hann');
        plot(w,h,'displayname',['N=' int2str(N(k))]);hold on
        title('Hanning');
end
legend
clear h f w
%%% All together    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
idx=1;
for k=1:length(N)
    
    [h(k,:),f,w(k,:)]=lowFltr(wc,fs,N(k),'hamm');
    eval(['f' num2str(idx) '= f;'])
    plot(w(k,:),h(k,:),'displayname',['Hamming N=' int2str(N(k))],'linewidth',2);hold on
    title('Hamming')
    idx=idx+1;
    [h(2+k,:),f,w(2+k,:)]=lowFltr(wc,fs,N(k),'hann');
    eval(['f' num2str(idx) '= f;'])
    plot(w(2+k,:),h(2+k,:),'--','displayname',['Hanning N=' int2str(N(k))],'linewidth',1.5);hold on
    
    idx=idx+1;
end
title('All filters');
legend
clear f k
%%% Filter on x(t)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%t=n*Ts
N=2^9;
Ts=1/fs;
n=0:(N-1);
%xt = sin(15*t) + 0.25*sin(200*t); 
% w1=15  f1=(15/2*pi)*Ts
% w2=200 f2=(200/2*pi)*Ts
% fs>2*max(f1,f2)
xn = sin(15*n*Ts) + 0.25*sin(200*n*Ts);
figure
[Yf,f]=defFft(fs,N,xn);
close;
names={'Hamming N=21' 'Hanning N=21' 'Hamming N=41' 'Hanning N=41'};
figure
for k=1:4
    subplot(2,2,k)
    
    eval(['filt=filter(f' int2str(k) ',1,xn);']);
%     Xf=fftshift(fft(filt));
    [Xf,faxis]=defFft(fs,length(filt),filt);hold off;
    plot(f,abs(Yf),'displayname','X(f)','linewidth',2);hold on;
    plot(f,abs(Xf),'r--','displayname','X(f) filtered','linewidth',2);
    title(names(k));
end

%%   fs=50   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wc=0.5*pi; fs = 50; N=[21 41];

%%% All together    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
idx=1;
for k=1:length(N)
    
    [h(k,:),f,w(k,:)]=lowFltr(wc,fs,N(k),'hamm');
    eval(['f' num2str(idx) '= f;'])
    idx=idx+1;
    [h(2+k,:),f,w(2+k,:)]=lowFltr(wc,fs,N(k),'hann');
    eval(['f' num2str(idx) '= f;'])
    idx=idx+1;
end
clear f k

N=2^9;
Ts=1/fs;
n=0:(N-1);
xn = sin(15*n*Ts) + 0.25*sin(200*n*Ts);
figure
[Yf,f]=defFft(fs,N,xn);
close;
names={'Hamming N=21' 'Hanning N=21' 'Hamming N=41' 'Hanning N=41'};
figure
for k=1:4
    subplot(2,2,k)
    eval(['filt=filter(f' int2str(k) ',1,xn);']);
    [Xf,faxis]=defFft(fs,length(filt),filt);hold off;
    plot(f,abs(Yf),'displayname','X(f)','linewidth',2);hold on;
    plot(f,abs(Xf),'r--','displayname','X(f) filtered','linewidth',2);
    title(names(k));
end
