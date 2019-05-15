function foo = g_square(t)
% g_square	square pulse function
%
% 1 where -0.5 < t <= 0.5, 0 otherwise

foo = (t<=0.5) & (t>-0.5);
