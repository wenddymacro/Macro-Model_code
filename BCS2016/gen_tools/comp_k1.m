function Gk1    = comp_k1(Gk,approx)

[np,na] = size(approx.coef);
np      = np-1;
cn      = approx.coef;
bnds    = approx.bnds;
Gk1     = zeros(size(Gk,1),na);
for i=1:na
    tmp     = 2*(log(Gk(:,i))-bnds(i,1))/(bnds(i,2)-bnds(i,1))-1;
    Gk1(:,i)= exp(cheb(tmp,np)*cn(:,i));
end
    