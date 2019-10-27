function [x,w,int]=gauss_herm(n,f,varargin);
%
% Compute the coefficients of Hermite polynomials using the recursion
%
% H(n+1)=2H(n)-2nH(n-1)
%
p0	= 1;
p1	= [2;0];
for i=1:n-1;
   p	= 2*[p1;0]-2*i*[0;0;p0];
   p0	= p1;
   p1	= p;
end
%
% Compute the roots
%
x	= sort(roots(p));
%
% Compute the weights imposing that integration is exact for lower order polynomials
%
A	= zeros(n);
A(1,:)	= ones(1,n);
A(2,:)	= 2*x';
for i=1:n-2;
   A(i+2,:)	= 2*x'.*A(i+1,:)-2*i*A(i,:);
end
w	= A\[sqrt(pi);zeros(n-1,1)];

if nargin>1;
   int=w'*feval(f,x,varargin{:});
end
