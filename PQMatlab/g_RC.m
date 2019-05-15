function gout = g_RC(t,beta,T,foo)
% g_RC: the raised-cosine pulse shaping function
%function gout = g_RC(t,beta,T,foo)
%
% t	time
% beta	roll-off factor
% T	symbol period (baud)
% 
% foo	optional flag. set to 1 to use exp(i*...) instead of cos.
%	Just to see what happens :->

% searles 2/2/98
% foo flag 16/7/98

if nargin==3, foo=0; end

if ~foo
  gout = sinc(t/T).*cos(pi*beta*t/T)./(1-4*beta*beta*(t/T).^2);
else
  % a stupid idea I thought I'd try:
  gout = sinc(t/T).*exp(i*pi*beta*t/T)./(1-4*beta*beta*(t/T).^2);
end
