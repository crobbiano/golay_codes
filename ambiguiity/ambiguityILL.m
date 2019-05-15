function Z = ambiguityILL(Ut1, Ut2, normalize, fftPoint)
% IMPORTANT: the logical 0 must be presented with -1

A=xcorril(Ut1, Ut2);
B=fft(A,fftPoint);
Z=abs(B);
clear A B;
p = size(Z);
lengthX = p(1);
halfLength=floor(lengthX(1)/2);
A1=Z(1:halfLength,:);
A2=Z(halfLength+1:lengthX,:);
X=cat(1,A2,A1);
Z=X;
clear A1 A2 X;

k = strcmpi(normalize,'true');
if k == 1 
    m = max (Z);
    mm = max (m);
    Z=Z./mm;
end

