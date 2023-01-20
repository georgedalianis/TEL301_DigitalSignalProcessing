function [convOut,L] = myConv(x,y)
%UNTITLED Summary of myConv
%   convOut = conv(x,y)
%   Detailed explanation goes here
L=1:length(x)+length(y)-1;

yWithPadding=[zeros(1,length(x)-1) flip(y) zeros(1,length(x)-1)];
convOut=zeros(0,length(yWithPadding));
for idx=1:(length(yWithPadding)-(length(x)-1))
    convOut(idx)=sum(x.*yWithPadding(idx:(idx+length(x)-1)));
end
convOut=flip(convOut)
% plot([1:L],convOut,'linewidth',2)
end

