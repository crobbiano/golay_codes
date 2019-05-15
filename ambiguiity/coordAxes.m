function [x,y] = coordAxes(freq, time, sampleRate, pulseLength)
% time = imput sequence Length (x axis)
% freq = the frequency axis length (y axis)
% sampleRate =  this is the chip length in a single pulse. if the pulse is [1 1 1 1] 
%            ----
%           |    |
%        ----    ---- 
%      then  sampleRate = 4. 
% 
%
% tb = pulse jength
% x axis  label = [tao/tb]
% -------------------------------------------------------------------
% ********  VERIFICATION  *****************
% CHECKED with:  
% barker Pulselength = 13, sampleRate =1;
% rectangular pulse Pulselength = 1, sampleRate =14;
% rectangular pulse Pulselength = 1, sampleRate =8;
% -------------------------------------------------------------------
% ********  REFERENCE FOR CHECKING  *****************
% COMPARED with "Radar Signals" Levanon N.,Eli Mozeson, Wiley 2004
% -------------------------------------------------------------------

half = floor((time)/2);
x = -half:half;
x = x/sampleRate;


half = (freq)/2;
y1=1:half-1;
y2 = (-half):0;
y = cat(2,y2,y1);
y = (y*sampleRate*pulseLength)/(freq); %normalize according to the pulse length

