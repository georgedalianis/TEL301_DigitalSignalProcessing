function [numZ,denZ]= MyButtFilter(passBand,stopband,delta_p,delta_s,fs,show)
%UNTITLED3 Summary of this function goes here
Wp = 2*pi*passBand; Ws = 2*pi*stopband;
%     [N,Wc] = buttord(Wp,Ws,delta_p,delta_s,'s');   % Gives minimum order of filter
%     [z,p,k] = butter(N,Wc);         % buttap(N);   Butterworth filter design
%     hs = zp2sos(z,p,k);            % Converts to second order sections
%     freqz(hs,1024,fs);
    %%
    
    [N,Wc] = buttord(Wp,Ws,delta_p,delta_s,'s');
    [z,p,k] = buttap(N);%   Butterworth filter rang N
    [n,d]=zp2tf(z,p,k);
    [num,den]= lp2lp (n,d,Wc); 
    f=linspace(0,fs/2,2048); % frequency space
    if show==1
     figure;
     freqs(num,den,fs)
     anal_freq=freqs(num,den,fs);
     title('Butterworth filter Contineous')
    end 
    % analog filter frequency 
    fil_anal=freqs(num,den,2*pi*f)
    [numZ,denZ]=bilinear(num,den,fs); % discrete system
    if show==1
     figure
     freqz(numZ,denZ)
     digit_freq=freqz(numZ,denZ,f,fs);
     title('Butterworth filter Discrete')
    
     figure
     plot(f,mag2db(abs(fil_anal)),'r.','linewidth',2.5,'Displayname','Contineous')
     hold on 
     plot(f,mag2db(abs(digit_freq)),'b--','linewidth',1,'Displayname','Discrete')
     title('Compare Butterworth filter Contineous and Discrete')
    end
end

