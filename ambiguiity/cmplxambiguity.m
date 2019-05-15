function Z = cmplxAmbiguity(normalize, fftPoint, sampleRate, Ut1, Ut2)
% THIS FUNCTION IS THE ENTRY POINT
% This function plots the monostatic ambiguity function of one or two signals (crossambiguiity function)
%
%
% Imput parameters:
% normalize = it normalises the function if it is 'true'
%
% fftPoint = it is the n-point DFT (for more see the function "fft" in the MATLAB functions reference )
%
%sampleRate =  this is the chip length in a single pulse (the sampling rate). 
%  if the signal is [0 1 0] four samples per pulse
%            ----
%           |    |
%        ----    ---- 
%      then  sampleRate = 4. 
% 
%  if the signal is [0 1 0] - one sample per pulse
%            -
%           | |
%          -   - 
%      then  sampleRate = 1. 
%
% Ut1 = the first input signal 
% Ut2 = the second input signal - if the signals are the same... only the first is acceptable
%
% -------------------------------------------------------------------
% ********  SAMPLE m-SYNTAX -  a barker code length 13  *****************
% -------------------------------------------------------------------
%
%  barker=[-1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1];
%  Z = cmplxambiguity('true', 256, 1, barker, barker);
%     . 'true'  -  normalize it
%     .  256    - 256 point of FFT
%     .  1      - the sample rate is 1 i.e. one bit one symbol
%     .  barker - the vector of the input signlal (
%                           if I want to find the ambiguiity function of one signal only the second imput signal could be skipped)
%                           if I want to find the cross ambiguiity function of two signals only the second imput signal should is obligated and have ot have same length as the previous)
%
%
%
% -------------------------------------------------------------------
% ********  VERIFICATION  *****************
% -------------------------------------------------------------------
% CHECKED with:  
% barker Pulselength = 13, sampleRate =1;
% rectangular pulse Pulselength = 1, sampleRate =14;
% rectangular pulse Pulselength = 1, sampleRate =8;
%
% -------------------------------------------------------------------
% ********  REFERENCE FOR CHECKING  *****************
% COMPARED with "Radar Signals" Levanon N.,Eli Mozeson, Wiley 2004
% -------------------------------------------------------------------
%
%
% author: Vladimir A. Kyovtorov, e-mail: vladimir_ak@yahoo.com, copyright 2006

if nargin < 4
    error('wrong number of arguments !!!');
end

if ~(strcmpi(normalize,'true')|strcmpi(normalize,'false'))
   normalize = 'false';
end
    
try
    Ut2;
catch
    Ut2 = Ut1;
end

Z = ambiguityILL(Ut1, Ut2, normalize, fftPoint);


% ********  PLOT the ambiguity function  *****************

% --------------------- axes estimating -----------------------------
ambigSize = size(Z);
pulseLength = length(Ut1);
[x,y] = coordAxes(fftPoint, ambigSize(2), sampleRate, pulseLength);

% --------------------- ploting -----------------------------
figure;
mesh(x, y, Z);
title('Complex ambiguity function');
xlabel('$\tau$/T', 'interpreter', 'latex');
ylabel('$\nu$T', 'interpreter', 'latex');
zlabel('|c($\tau$,$\nu$)|', 'interpreter', 'latex')



