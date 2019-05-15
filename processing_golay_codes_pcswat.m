%% read files and match
clear all; stave = convertPCSWAT('C:\Users\iamchris\Documents\ONR_files\PQMatlab\staves.txt');
plot(20*abs(stave.elements.data))
incident = convertIncident('C:\Users\iamchris\Documents\PC SWAT 10.5.1\Win64\SignalIncident.txt');
incident = incident(:,2);
matched = matchToIncident(stave.elements.data, incident);
figure(1); subplot(3,1,1);plot(20*abs(stave.elements.data));subplot(3,1,2); plot(abs(incident)); subplot(3,1,3); plot(abs(matched))

%% look at fft
Fs = 100000;            % Sampling frequency                    
T = 1/Fs;               % Sampling period   
c = 1500;               % sound speed

X = stave.elements.data;    
L = length(X);             % Length of signal
t = (0:L-1)*T;        % Time vector
Y = fft(X);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(2)
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
%% add in doppler
theta = linspace(-pi,pi,3143);
thradPRI = theta;
% M = zeros(length(X),length(thradPRI));
M = matched(:)*exp(1i*thradPRI);

nmf_1 = 1;
nmf_2 = length(X);
minM = 3e-6;
figure(3)
imagesc(theta, (nmf_1:nmf_2)/Fs*c/2,20*safelog(abs(M(nmf_1:nmf_2,:))/max(max(abs(M))),10,10*minM))
set(gca,'YDir','normal')
xlabel('Doppler (rad)');
ylabel('Range (m)');
colormap(jet), colorbar;