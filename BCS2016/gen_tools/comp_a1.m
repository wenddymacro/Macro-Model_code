function Ga1    = comp_a1(Ga,abar,approx)

[np,nz] = size(approx.cofn);
np      = np-1;
cn      = approx.cofn;
cc      = approx.cofc;
bnds    = approx.bnds;
Ga1     = zeros(size(Ga,1),nz);
for i=1:nz
    I       = find(Ga(:,i)<abar(i));
    J       = find(Ga(:,i)>=abar(i));
    tmp     = 2*(log(Ga(I,i))-bnds(i,1))/(bnds(i,2)-bnds(i,1))-1;
    Ga1(I,i)= exp(cheb(tmp,np)*cn(:,i));
    tmp     = 2*(log(Ga(J,i))-bnds(i,2))/(bnds(i,3)-bnds(i,2))-1;
    Ga1(J,i)= exp(cheb(tmp,np)*cc(:,i));
end
    