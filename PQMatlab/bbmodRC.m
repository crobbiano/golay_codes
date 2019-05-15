% basebandmodRC     baseband modulation with raised cosine PSF
%
%[p,ts] = basebandmodRC(code1,smp,rolloff)
%
% code1 - code sequence
% smp   - upsampling factor
% roff  - rolloff for raised cosine, may omit.
%
% p     - BB modulated code
% t     - time vector (in samples) of bb signal,
%         0 at sample point of first element.

% SJS 8/3/07

function [p,ts] = basebandmodRC(code1,smp,rolloff)

if nargin<3 || isempty(rolloff)
    rolloff = 0.9;
end

N = length(code1);

%ts = -smp:(N)*smp;  % was N+1
ts = -2*smp:(N+1)*smp;
p = zeros(size(ts));
for q = 1:N
    grc = g_RC(ts-(q-1)*smp,rolloff,smp);
    p = p + code1(q)*grc;
end
