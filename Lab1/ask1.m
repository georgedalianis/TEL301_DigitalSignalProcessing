clc; close all; clear all;
%%  1.A Convolution
n=15
x=round(10*rand(1,n));
y=round(10*rand(1,n));
figure('name'," Signals")
plot([1:n],x,'displayname',"x(n)");hold on
plot([1:n],y,'displayname',"y(n)");hold on
title("Two random signals x(n) and y(n) ");

figure('name'," Compare convolutions")
subplot(3,1,1)
[myconv,n] = defConv(x,y);
plot(n,myconv);
title("Convolution by definition [arrays]")
subplot(3,1,2)
convXY=conv(x,y);
plot(n,convXY);
title("Convolutionwith  with Matlab conv()")
subplot(3,1,3)
[out,l]= myConv(x,y);
plot(l,out);
title("Convolution with our funtion myConv()")

%% 1.B Fourier
% clc;close all;
X_f = fft(x,length(l)); % X(f)=fourier tranform of x(t)
Y_f = fft(y,length(l));
con=ifft(X_f.*Y_f);
figure('Name'," Convolution in time and in frequency ")
subplot(2,1,1)
plot(n,con); % Show  matlab convolution
title("Convolution with F^{-1}[X_f(f)*Y_f(f)] ")
subplot(2,1,2)
plot(n,convXY);
title("Convolution in time with Matlab conv()")
figure
plot(n,con,'displayname',"fft",'linewidth',5);
hold on
plot(n,convXY,'displayname',"conv",'linewidth',3);
title("Convolutions comparison fft vs conv ")

%% 2
%clc;close all;clear all;
figure('Name',"")
t=[0:0.001:0.5];    % time
f1=12;
f2=1.5/2;
fs=2*max(f1,f2);                     % Nyquist

x=5*cos(24*pi*t)-2*sin(1.5*pi*t);    % signal

figure(10);
plot(t,x);

% Sampling Nyquist
t1=(0:(1/48):0.5);
subplot(2,2,1)
y48=showSig(t,x,t1,5*cos(2*f1*pi*t1)-2*sin(2*f2*pi*t1),"T_s=1/48s");
% Sampling T_s=1/24ss
%fs=f1;
t1=(0:(1/24):0.5);
subplot(2,2,2)
y24=showSig(t,x,t1,5*cos(2*f1*pi*t1)-2*sin(2*f2*pi*t1),"T_s=1/24s");
% Sampling T_s=1/12s
fs=12
t1=(0:(1/fs):0.5);
subplot(2,2,3)
y12=showSig(t,x,t1,5*cos(2*f1*pi*t1)-2*sin(2*f2*pi*t1),"T_s=1/12s");
% Sampling T_s=1/A s
fs=47
t1=(0:(1/fs):0.5);
subplot(2,2,4)
y47=showSig(t,x,t1,5*cos(2*f1*pi*t1)-2*sin(2*f2*pi*t1),"T_s=1/A s");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('name',"Fourieres")
subplot(2,2,1)
N=length(y48);
[outfft,faxis] =defFft(48,N,y48);
subplot(2,2,2)
N=length(y24);
[outfft,faxis] =defFft(24,N,y24);
subplot(2,2,3)
N=length(y12);
[outfft,faxis] =defFft(12,N,y12);
subplot(2,2,4)
N=length(y47);
[outfft,faxis] =defFft(47,N,y47);

%% 3.A
%clc;close all; clear all;
figure
t=[0:0.001:0.5];    % time
x=10*cos(2*pi*20*t)-4*sin(2*pi*40*t+5);
title("3.A signal")
fs=2*max(20,40);
% N=128 Samplings fs=255
N=128;
fs=255; % Ts=(0.5-0)/128
ty=(0:(1/fs):0.5);
xt=showSig(t,x,ty,(10*cos(2*pi*20*ty)-4*sin(2*pi*40*ty+5)),"3.A with T_s=1/255_s");
figure('name',"3.A Frequency Specrum")
[outfft,faxis] =defFft(fs,N,xt);

%% 3.B 100:125:475
%clc;close all; clear all;
t=[0:0.0001:0.3];    % time
f0=(100:125:475);   % signal frequencies
phi=47;             % phi = team number
fs=8000;            % sampling fryquency
figure
pos=1;    % position in figure

for i=1:length(f0) % for each signal frequency
    x=sin(2*pi*f0(i)*t+phi); % calculate signal
    subplot(length(f0),2,pos)
    plot(t,x);title("Signal x(t) with f_0="+f0(i)) % show signal
    subplot(length(f0),2,pos+1)
    [outfft,faxis] =defFft(fs,length(t),x);title("Fourier X(f) f_0="+f0(i)); % show Fourier
    tbf0(i,:)=outfft;
    foAxis(i,:)=faxis;
    pos=pos+2;
end
% sampling
x=0.*t;
fs=8000;            %  n=t*Ts => t=n/fs
figure
for sg=1:length(f0) % for each signal frequency
    for i=1:length(t)
        x(i)=sin(2*pi*(f0(sg)/fs)*i+phi); % x(t)=sin(2*pi*f0*(n/fs)+phi)=sin(2*pi*(f0/fs)*n+phi)
    end
    subplot(length(f0),2,(sg*2)-1)
    plot(t,x);
    title("signal with f_0="+f0(sg)+" and "+fs)
    subplot(length(f0),2,sg*2)
    [outfft,faxis] =defFft(fs,length(t),x);
    
    tbfs(sg,:)=outfft;
    fsAxis(sg,:)=faxis;
end
figure('name',"compare fouriers")
pos=1;
for row=1:length(f0)
    subplot(length(f0),2,pos);pos=pos+1;
    plot(foAxis(row,:),tbf0(row,:))
    title("f_0="+f0(row))
    subplot(length(f0),2,pos);pos=pos+1;
    plot(fsAxis(row,:),tbfs(row,:))
    title("f_0="+f0(row)+"f_s=8000")
end

%% 3.B  7525:125:7900

%clc;close all; clear all;
t=[0:0.0001:0.3];    % time
f0=(7525:125:7900);   % signal frequencies
phi=3*47;             % phi = team number
fs=8000;            % sampling fryquency
figure
pos=1;    % position in figure

for i=1:length(f0) % for each signal frequency
    x=sin(2*pi*f0(i)*t+phi); % calculate signal
    subplot(length(f0),2,pos)
    plot(t,x);title("Signal x(t) with f_0="+f0(i)) % show signal
    subplot(length(f0),2,pos+1)
    [outfft,faxis] =defFft(fs,length(t),x);title("Fourier X(f) f_0="+f0(i)) % show Fourier
    pos=pos+2;
    tbf0(i,:)=outfft;
    foAxis(i,:)=faxis;
end
% sampling
x=0.*t;
fs=8000;            %  n=t*Ts => t=n/fs
figure
for sg=1:length(f0) % for each signal frequency
    for i=1:length(t)
        x(i)=sin(2*pi*(f0(sg)/fs)*i+phi); % x(t)=sin(2*pi*f0*(n/fs)+phi)=sin(2*pi*(f0/fs)*n+phi)
    end
    subplot(length(f0),2,(sg*2)-1)
    plot(t,x);
    title("signal with f_0="+f0(sg)+" and "+fs)
    subplot(length(f0),2,sg*2)
    [outfft,faxis]=defFft(fs,length(t),x);
    tbfs(sg,:)=outfft;
    fsAxis(sg,:)=faxis;
end

figure('name',"compare fouriers")
pos=1;
for row=1:length(f0)
    subplot(length(f0),2,pos);pos=pos+1;
    plot(foAxis(row,:),tbf0(row,:))
    title("f_0="+f0(row))
    subplot(length(f0),2,pos);pos=pos+1;
    plot(fsAxis(row,:),tbfs(row,:))
    title("f_0="+f0(row)+"f_s=8000")
end
