function c = xcorril(a,b)

%if nargin == 1
%	b = a;
%end

[ma,na] = size(a);
[mb,nb] = size(b);

b = conj(b(mb:-1:1,:));
apad = [a ; zeros(mb-1,na)];

c = zeros(na+nb-1,na+nb-1);

for k=1:(na+nb-1)
	count = k		*(k<min(na,nb)) ...
		+min(na,nb)	*(k>=min(na,nb))*(k<=max(na,nb)) ...
		+(na+nb-k)	*(k>max(na,nb));

	starta = 1		*(k<=nb) ...
		+(k-nb+1)	*(k>nb);

	startb = (nb-k+1)	*(k<=nb) ...
		+1		*(k>nb);

	for i=0:(count-1)
		c(i+1,k) =  filter( b(:,startb+i), 2, apad(:,starta+i) ); 
	end
end

