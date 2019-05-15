function gH = golay_Hmat(arg,lgc)
%

if isscalar(arg)
    N = arg;
    ovec = ones(1,N-1);
else
    ovec = arg;
    N = 2^(length(ovec));
end
if nargin<2 | isempty(lgc)
    lgc = false;
end

l2N = log2(N);
if l2N~=round(l2N)
    error('input arg should be power of 2')
end

H = logical(1);

for q = 1:l2N
    
    %Htilde = circshift(H,size(H,1)/2);  % 2^(q-2) ?
    Htilde = circshift(H,floor(2^(q-2)));  % 2^(q-2) ?
    % circshift of 1 is same as 0 when 
    if ovec(q)
        Htilde = ~Htilde;
    end
    H = [H Htilde;
         H ~Htilde ];
     
end

if lgc
    gH = H;
else
    gH = 2*H-1;
end
