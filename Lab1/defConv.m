function [convOut,L] = defConv(x,y)
%UNTITLED2 Summary of defConv
%
%   Detailed explanation goes here
%%
N=length(x)+length(y)-1;
L=1:N;
convOut=0*zeros(1,N); 
for i=1:N
    for k=1:length(x)
        if(i-k+1>0 & i - k + 1 <= length(y))
            convOut(i)= convOut(i) +(x(k).*y(i-k+1));
        end
    end
end
end

