% safelog   safeguarded logarithm
%
% prevent logs of zero/negative values
%
% dat   data matrix
% base  10, 2, or otherwise for base e
% tol   minimum value.  If omitted, will use smallest nonzero datum.

function lout = safelog(dat,base,tol)

if nargin<2 || isempty(base)
    base = 'e';
end
if nargin<3 || isempty(tol)
    dval = sort(dat(:));
    foo = find(dval>0);
    tol = dval(foo(1));
end

duse = max(dat,tol);

switch base
    case 2
        lout = log2(duse);
    case 10
        lout = log10(duse);
    otherwise
        lout = log(duse);
end
