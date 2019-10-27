function Gk1    = comp_kp(Gk,approx,i)

np      = size(approx.coef,1);
np      = np-1;
cn      = approx.coef;
bnds    = approx.bnds;
tmp     = 2*(log(Gk)-bnds(i,1))/(bnds(i,2)-bnds(i,1))-1;
Gk1     = exp(cheb(tmp,np)*cn(:,i));
    