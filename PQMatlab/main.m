% main. m      
% Emulation of radar images of point targets in scene. 
% Use binary sequence to schedule the transm. of Golay compl. waveforms.
% Use real sequence to weight returns in each PRI. 
% *** version 3 ***
%
% required coding of sim_code_clut_PRI.m

% W Dang July 2010

disp('   *** version 3 ***')

clear all

addpath([pwd, '\Utility'])
%% system params

c = 299792458;
if 0
    fc = 1e9;
    %fo = 50e6;
    omgc = 2*pi*fc;
    T = 1e-2;  % PRI
    gnch = 8;
    sampF = 10e6;
    sampT = 1/sampF;
    ups = 8;
    
    ncode = 2;
    npow = 2;
    desN = [];
elseif 0
    fc = 17e9;        % carrier frequency
    omgc = 2*pi*fc;
    prf = 20e3;       % pulse repetition frequency
    T = 1/prf;        % PRI

    sampF = 50e6;
    sampT = 1/sampF;

    gnch = 64;        % # chips in Golay sequence
    Tch = 100e-9;     % time of one chip
    Tg = gnch*Tch;    % time of whole Golay code
    nsg = Tg*sampF;   % # samples required for Golay code
    ups = nsg/gnch;   % upsampling factor
    % may need to round the above, thereby changing values above

    ncode = 2;        % Golay pair
    desN = 16;        % pulse train length 
    
%     seqmethod = 4;    % select combination of (P,Q)
    seqmethod = 3;    % select combination of (P,Q)
else 
    % sonar
    c = 1500;
    fc = 100e3;        % carrier frequency
    omgc = 2*pi*fc;
    prf = 40;       % pulse repetition frequency
    T = 1/prf;        % PRI

    sampF = 100e3;
    sampT = 1/sampF;

    gnch = 64;        % # chips in Golay sequence
    Tch = 5e-5;     % time of one chip
    Tg = gnch*Tch;    % time of whole Golay code
    nsg = Tg*sampF;   % # samples required for Golay code
    ups = nsg/gnch;   % upsampling factor
    % may need to round the above, thereby changing values above

    ncode = 2;        % Golay pair
    desN = 16;        % pulse train length 
    
%     seqmethod = 1;    % select combination of (P,Q)
    seqmethod = 3;    % select combination of (P,Q)
end

%% targetparams

RPS = sampT*c/2;  % for info... range per sample
switch 5
    case 1
        scat = [
            1       2/c*1000    0
%             1       2/c*1010    0
%             1       2/c*1020    0
%             1       2/c*1030    0
%             1       2/c*1040    0
%             1       2/c*1050    0
%             1       2/c*1060    0
%             1       2/c*1070    0
            .03       2/c*1110    1.5 %change from 0.033
            1       2/c*1500    0
            .03       2/c*1600    -1.5 %change from 0.033
            1       2/c*1700    0
            ];
    case 2
        scat = [1 0 10];% 1 2/c*20*RPS 0];
    case 3
        scat = [1 2/c*1000 0];
    case 4
        % sonar 1
        scat = [
            1       2/c*100    1.5
            .03       2/c*110    -1
            1       2/c*120    0
            .03       2/c*130    0
            ];
    case 5
        % sonar 2
        scat = [
            .8       2/c*98    0
            1       2/c*100    1.5
            .3       2/c*100.5    0
            .03       2/c*101    -1
            ];
end
nscat = size(scat,1);
if 1
    % random phases on returns
    scat(:,1) = scat(:,1).*exp(i*2*pi*rand(nscat,1));
else
    disp('no random phases on returns')
end


% frequency shift range arguments
vmax = max(abs(scat(:,3)));
farg = [sampF/10000/10 50*vmax/c*fc];  
disp('frequency arguments:'),  disp(farg)


%% Transmit waveform and receive filter

% Golay pair 
code = golay_Hmat(max(gnch,ncode)); code = code(1:ncode,:);

% sequence of code indices
seq  = zeros(1,desN);
qseq = zeros(1,desN);
switch seqmethod
    case 1
        % alternating
        designmet = 'alternating';
        seq = mod(0:desN-1,ncode);
        nullord = 0; 
        qseq = ones(1,desN); 
    case 2
        % PTM
        designmet = 'PTM';
        nullord = log2(desN) - 1; 
        seq = tmcomp(nullord,ncode);
        qseq = ones(1,desN); 
    case 3
        designmet = 'binomial';
        seq = mod(0:desN-1,ncode); 
        nullord = desN-2; 
        for priN = 1 : desN
            qseq(priN) = nchoosek(desN-1,priN-1); 
        end
    case 4
        % "shuffled" (random but equal proportions)
        designmet = 'shuffled';
        seq = mod(randperm(desN)-1,ncode);
        nullord = 0; 
        qseq = ones(1,desN);
end
disp(['  ', designmet]); 
Norg = desN;

% Baseband radar return and matched filter
bmth = 'r';

for q = 1:ncode
    ix = find(seq==q-1);
    dv = (ix-1)*T;
    [bd,b] = sim_code_clut_PRI(code(q,:),ups,bmth,scat,fc,dv,sampF);
    if q==1
        ret = zeros(desN,length(bd));
        bmat = zeros(desN,length(b));
    end
    ret(ix,:) = bd;  % baseband return
    bmat(ix,:) = b(ones(desN/ncode,1),:); % baseband matched filter
end


%% receive processing (zero Doppler)

nmf = size(ret,2)+size(bmat,2)-1;
mf = zeros(desN,nmf);

theta = linspace(-pi,pi,3143);
thradPRI = theta;
M = zeros(nmf,length(thradPRI));

% This takes the matched filter output and shifts each sample (delay) by a 
% phase (in theta variable) then sums them all together over the delays in
% a single code (the desN)
for q = 1 : desN
    mfout = qseq(q)*conv(ret(q,:),conj(bmat(seq(q)+1,end:-1:1)));
    mf(q,:) = mfout;
    Mq = mfout(:)*exp(sqrt(-1)*(q-1)*thradPRI);
    M = M + Mq;
end

% observation window
flag = 1;
row_index = 0;
while flag
    row_index = row_index + 1;
    if max(abs(M(row_index,:))) > 0
        flag = 0;
    end
end
nmf_1 = row_index;
flag = 1;
row_index = nmf+1;
while flag
    row_index = row_index - 1;
    if max(abs(M(row_index,:))) > 0
        flag = 0;
    end
end
nmf_2 = row_index;

% Plotting radar image
minM = 3e-6;
figure;
% imagesc(theta, (nmf_1:nmf_2)/sampF,20*safelog(abs(M(nmf_1:nmf_2,:))/max(max(abs(M))),10,10*minM))
% ylabel('Delay (sec)');
imagesc(theta, (nmf_1:nmf_2)/sampF*c/2 - 2.5,20*safelog(abs(M(nmf_1:nmf_2,:))/max(max(abs(M))),10,10*minM))
set(gca,'YDir','normal')
xlabel('Doppler (rad)');
ylabel('Range (m)');
colormap(jet), colorbar;
title({['N = ', num2str(desN), ' pulses, order of null: M = ', num2str(nullord), ', ' ]; [designmet, ' design, ', 'output in dB']});

%% End 