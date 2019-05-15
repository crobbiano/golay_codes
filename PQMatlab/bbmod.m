% basebandmod     take one or 2 codes and form baseband PCM signals.
%
% code1     a code.  
%               binary:  logical -OR- [-1 1] valued.
%               other:  complex phasors for PCM eg [1 i -1 -i]
% smp       upsampling factor.  integer.
% filt      [cutoff order]  as per fir1.  Or -1 for no filter.
% dotrim    true to trim filter residuals, false to leave. (OPTIONAL)
%
% p     the code(s) upsampled and filtered prior to modulation

% SJS 8/6/06
% based heavily on Sofia's code
%
% 15/8/06
% further converted from code2baseband.m, from the multiplexed Golay code
% experiments.
% 
% 25/8/06
% can now supply optional arg to include/exclude filter residuals
%
% 9/3/07
% fixed bug:  dotrim == true was NOT trimming and vice versa
% also added n2=floor(n/2) in trim bit.
%
% 30/3/07
% no filtering done if filter argument is -1

function [p] = basebandmod(code1,smp,filt,varargin)
% former parameters excised:
% code2     as above or []
% coding    'IQ', 'freq', 'time'
% also output parameter w removed.
% w     complex baseband modulated sequence
code2 = [];
% leaving code2 stuff /in situ/ in case I need it later

% default parameters
if isempty(filt)
    filt = [.25 16];
end
if isempty(smp)
    smp = 1;
end
lva = length(varargin);
if lva>0
    dotrim = logical(varargin{1});
else
    dotrim = false;
end

% arrange codes
clen = length(code1);
codes = [code1(:) code2(:)]';  % each row is a code
if islogical(codes) || (~any(any(codes<0)) && max(max(imag(codes)))==0)
    % second || clause was:  any(any(codes==0))
    % however we might want 0 values if code is already phase-coded
    % so instead check that code is real-valued and has no negative values
    codes = 2*double(codes)-1;
end

% upsample codes
codes = kron(codes,ones(1,smp));

% filter upsampled codes
if filt(1)~=-1 % 30/3/07
    
    c = filt(1);
    if isscalar(filt)
        n = 16;
    else
        n = filt(2);
    end
    %if license('checkout','signal_toolbox')
    %    b=fir1(n,c,blackmanharris(n+1));
    %else
    b=fir1(n,c);  % using jolly rogered version
    %end
    a=1;
    [p,zf] = filter(b,a,codes');
    %!@#$
    % pick one or other of these based on subjective desire
    if ~dotrim  % ~ prepended 9/3/07
        %if lva==0,
        disp('including filter residuals');
        %end
        p = [p; zf]';
    else
        %if lva==0,
        disp('trimming filter residuals');
        %end
        n2 = floor(n/2);  % 9/3/07
        p = [p(n2+1:end,:); zf(1:n2,:)]';
    end

else
    % filter suppressed
    p = codes;
    disp('basebandmod: no filter')
end
