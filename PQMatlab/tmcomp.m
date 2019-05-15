% tmcomp    Generate Prouhet-Thue-Morse sequence for arbitrary base
%
% Generate t, the sequence of integers on (1, base) 
% such that the sum of all (powers of) integers i (mod base)
% on (0,N-1) (mod base)
% having same t(i)
% is the same for all possible t.
%
% bpow - maximum power of integers for which sum-property holds
% base - as above.  defaults to 2 (binary case: original PTM seq)

% SJS 30/3/07

function t = tmcomp(bpow,b)

if nargin<2 | isempty(b)
    b = 2;
end

N = b^(bpow+1);
%fprintf('base %i,  N %i\n',b,N)

% this works for b==2
num = (0:N-1);  % numbers from 0 to N-1
s = zeros(1,N); % number of 1s in binary representation
                % "binary digit sum" (for b==2 only)
while any(num>=1)
    m2 = mod(num,b);
    s = s + m2;
    num = floor(num/b);
end
t = mod(s,b);    % parity (odd number of 1s) if b==2
% t is (generalised?) P-T-M seq?

% could've generated t via production rules
% 0->01  1->10
% what about higher bases?
% 0->01  1->23  2->12  3->30
%
% or via concatenation...???  higher bases?

