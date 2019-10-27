function approx=smoothdecrule(a,a1,abar,amin,amax,np)
%
% approx=smoothdecrule(a,a1,abar,amin,amax,np)
%   Input:  a   : vector of initial assets a(t)
%           a1  : vector of nex period assets a(t+1)
%           abar: Threshold values
%           amin: minimal value for assets in the grid
%           amax: maximal value for assets in the grid
%           np  : order of polynomial approximation
%
%   Output: approx structure
%           approx.cofn : Coefficients of the approximation in normal times
%           approx.cofc : Coefficients of the approximation in crisis times
%           approx.bnds : bounds for transformation function mapping asset
%                         values into [-1,1]
%
nz  = size(a1,2);
cn  = zeros(np+1,nz);
cc  = zeros(np+1,nz);
b   = zeros(nz,3);
an  = log(amin);
bc  = log(amax);
for j=1:nz
    I       = find(a(:,j)<abar(j));
    J       = find(a(:,j)>=abar(j));
    bn      = log(abar(j));
    ac      = log(abar(j));
    x0      = 2*(log(a(I,j))-an)/(bn-an)-1;
    x1      = 2*(log(a(J,j))-ac)/(bc-ac)-1;
    c0      = cheb(x0,np)\log(a1(I,j));
    c1      = cheb(x1,np)\log(a1(J,j));
    
    cn(:,j) = c0(:);
    cc(:,j) = c1(:);
    b(j,:)  = [an bn bc];
end

approx.cofn= cn;
approx.cofc= cc;
approx.bnds = b;

