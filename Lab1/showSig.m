function sig = showSig(tx,x,ty,y,ttl)
%% UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%%
sig=y;
plot(tx,x,'linewidth',1,'displayName',"signal");hold on
stem(ty,y,'linewidth',1,'displayName',"samples");hold on
title(ttl);
legend
end

