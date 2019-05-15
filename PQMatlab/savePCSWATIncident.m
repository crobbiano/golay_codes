
for ii=1:321
    sig(ii) = b(ii+7);
end

%% write to file

times = 0:1e-5:0.0032;
reall = sig;
imagg = zeros(size(sig));

%% 
fileID = fopen('signalPCSWAT.txt','w');
for k=1:length(sig)
    fprintf(fileID,'%f %f %f\n',times(k), reall(k), imagg(k));
end
fclose(fileID)