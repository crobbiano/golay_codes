% sim_code_clut_delay
% function [ret,bb] = sim_code_clut_PRI(g1,smp,method,scat,fc,dvec,sampf)
% simulate radar returns (at baseband) from multiple scatterers
% with Golay-coded signal.
% actually could be any code, need not be Golay.
% phase is tweaked to account for delay-dependent Doppler phase
% at different pulse intervals.
%
% gpair     1xN     (Golay) code
% smp       scalar  upsampling factor
% method    scalar (for now)
%                   'f'     upsample & lowpass filter
%                   'r'     use raised cosing pulse
% scat      Kx3     scatterer information
%                   each row is a 3-tuple:
%                   RCS     complex reflection coeff
%                   delay   round trip delay (within PRI)
%                   speed   velocity 
% fc        scalar  carrier frequency, same for each pulse
%           1xV     carrier frequency, one for each pulse
% dvec      1xV     vector of pulse delays.  typically (0:foo-1)*PRI
%                   these are delays of initial sample:  signals themselves
%                   not translated, just have aptly modified phase.
% sampf     scalar  sampling frequency
%
% ret       VxL     returned & demodulated signals, one per row.
%                   one row for each delay.
% bb        1xL     baseband signals, pre-transmission
%                   first & second codes respectively.

% SJS 9/3/07
% SJS 13/3/07   second code now synthesised from first if not supplied.
%               should really suppress second-code output instead.
% SJS 22/3/07   time-delay of second code
%
% SJS 3/4/07    ripped off from sim_coded_cluttered for use with
%               Prouhet-Thue-Morse ordering of Golay pair
% SJS 25/9/07   modification allowing fc to be a vector
%               so frequency can be wobbled across pulses

function [ret,bb] = ...
    sim_code_clut_PRI(g1,smp,method,scat,fc,dvec,sampf)

npi = length(dvec);

if isscalar(fc)
    fc = fc(1,ones(1,npi));
elseif length(fc)~=npi
    error('fc must be scalar, or vector of same length as dvec.')
end
omgc = 2*pi*fc;
sampT = 1/sampf;
c = 299792458;
% c = 1500; % FIXME

nscat = size(scat,1);
maxdel = max(scat(:,2));

% zero-delay baseband signals
switch method(1)
    case 'f'
        fparm = [0.5 smp];
        dotrim = false;
        b1 = bbmod(g1,smp,fparm,dotrim);  % basebandmod was renamed
    case 'r'
        roll = 2.1;%0.9;
        b1 = bbmodRC(g1,smp,roll);  % basebandmodRC was renamed
    case 'u'
        % added 5/9/07, needs checking:
        b1 = bbmod(g1,smp,-1);
    otherwise
        error('no such method: %s',method);
end
bb = b1;

nb = length(b1);  % should be same as b2
maxtn = 100*ceil((maxdel*sampf+nb)/100);
ret1 = zeros(nscat,maxtn);
maxt = (maxtn-1)*sampT;  % ?
td = (0:maxtn-1)*sampT;

% compose baseband signals by cumulatively summing delayed versions
% (one for each scatterer)
for q = 1:nscat
    a = scat(q,1);
    d = scat(q,2);
    v = scat(q,3);  
    
    d1 = delay_by_tau(sampf,b1,d,maxt);
    
% since fc may now vary from pulse to pulse, must deal with phase entirely
% in the pulse-loop below.
%    if foobar
%        ret1(q,:) = a*d1 * exp(-i*omgc*d) ...
%            .* exp(i*(omgc*2*v/c)*(0:sampT:maxt));
%    else
        ret1(q,:) = a*d1;
%    end
end

% account for PRI Doppler phase
c = 299792458;
% c = 1500; % FIXME
ret = zeros(npi,maxtn);
vel = scat(:,3);
del = scat(:,2);
for q = 1:npi

%    if foobar
%        % account only for interpulse phase effects
%        fayze = 2*vel/c*dvec(q);
%        foo = repmat(exp(i*omgc*fayze),[1 maxtn]) .* ret1;
%    else
        % account for interpulse, and transmit-receive delay/Doppler
        trm = 2*vel/c*(dvec(q)+td);
        fayze = trm-del(:,ones(1,maxtn));
        foo = exp(i*omgc(q)*fayze) .* ret1;
%    end
    
    %fayz1(q) = fayze(1);  % for debuggerising
    
    ret(q,:) = sum(foo,1);

end

%ret=ret+.01*randn(size(ret,1),size(ret,2));

%disp(fayz1)

return
