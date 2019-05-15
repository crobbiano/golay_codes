% delay_by_tau  delay a pulse by tau secs.  may be non-integral #samples
%
% signal is zero outside of domain of original pulse.
%
% targ      vector of sample times OR sampling frequency
% p         pulse to be delayed.
% tau       # secs to delay
% maxt      maximum sample time of return vector.  May omit.
%           but best not to omit!
%
% pdel      delayed pulse
% tbig      sample times of delayed pulse

% SJS 16/6/06
% SJS 9/3/07 hacked to make initial sampling time = 0.
%            hopefully this won't screw with existing code.
%            places of change are commented with date.

function [pdel,tbig] = delay_by_tau(targ,p,tau,maxt)

if isscalar(targ)
    sampf = targ;
    Nsamp = length(p);
    %tvec = (1:Nsamp)/sampf;
    % SJS 9/3/07:
    tvec = (0:Nsamp-1)/sampf;
else
    tvec = targ;
    sampT = tvec(2)-tvec(1);
    sampf = 1/sampT;
    Nsamp = length(tvec);
end

if nargin<4 || isempty(maxt)
    %mingap = 0.1;
    %maxt = ceil(max(tvec)+tau+mingap);
    %maxt = max(tvec)+max(tau,0)+0.1*Nsamp/sampf;
    %maxt = Nsamp/sampf;  % i.e. limit to current size.
                         % too bad if that's not big enough
    % SJS 9/3/07:
    maxt = (Nsamp-1)/sampf;
end
%if maxt < max(tvec)+tau
%    error('desired maximum t is too small to encompass delayed pulse')
%end

%tbig = (1:maxt*sampf)/sampf;
% SJS 9/3/07:
tbig = (0:round(maxt*sampf))/sampf;
% note round needed 'coz of numerical errors.
% 0:1/sampf/maxt worked just fine!

imeth = 'spline';
%fprintf('delay_by_tau using interp1 with method %s\n',imeth);
pdel = interp1(tvec,p,tbig-tau,imeth,0);
