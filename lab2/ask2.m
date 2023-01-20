clc;
clear all;
close all;
%% 1 B
syms z 
NUM=[0.2 0];
DEN=[1 -0.7 -0.18];
display('------  1.B. H(f) z transform   ------')
syms z
H=tf(NUM,DEN);
sys = tf([0 0.2],[1 -0.7 -0.18],0.1,'Variable','z^-1')
zplane(NUM,DEN);
title('1.B zeros-poles')
%% 1 D
figure
f=[(-pi):(pi/128):(pi)];
freqz(NUM,DEN,f)
title('Response with frequency space');
figure
freqz(NUM,DEN)
title('Response without frequency space');

%%  1 E add pole z=1
NUM=[0 0.2 0];
DEN=[1 -1.7 0.52 0.18]; % s^3- 1.7s^2 0.52s + 0.18
display('------  1.E H(f) z transform   ------')
H=tf(NUM,DEN,0.1,'Variable','z^-1')
[R,P,K] = residuez(NUM,DEN);
figure
zplane(NUM,DEN);
figure
f=[(-pi):(pi/128):(pi)];
freqz(NUM,DEN,f)
title('Response with extra pole z=1');

%%  2 A
syms z n
display('------   H(f) z transform   ------')
pretty((4-3.5*z^(-1))/(1 -2.5*(z^(-1)) + z^(-2) ))
[R,P,K] = residuez([4 -3.5],[1 -2.5  1] );
display('------   Simple fractions   ------')
pretty(((4-3.5*z^(-1))/(1 -2.5*(z^(-1)) + z^(-2) ))==(R(1)/(1-(P(1))*z^-1))+(R(2)/(1-(P(2))*z^-1)))

%% 2 B
display('------   Invers z transform   ------')
x_n=iztrans((R(1)/(1-(P(1))*z^-1))+(R(2)/(1-(P(2))*z^-1)))
