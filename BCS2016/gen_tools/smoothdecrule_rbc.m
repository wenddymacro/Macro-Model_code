function approx=smoothdecrule_rbc(k,k1,kmin,kmax,np)
[~,na] = size(k1);
cn      = zeros(np+1,na);
b       = zeros(na,2);
for j=1:na
    a0      = log(kmin);
    b0      = log(kmax);
    x0      = 2*(log(k(:,j))-a0)/(b0-a0)-1;
    c0      = cheb(x0,np)\log(k1(:,j));
    
    cn(:,j) = c0(:);
    b(j,:)  = [a0 b0];
end

approx.coef= cn;
approx.bnds = b;

