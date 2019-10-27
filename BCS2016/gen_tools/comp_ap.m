function Ga1    = comp_ap(Ga,abar,approx,iz)

np      = size(approx.cofn,1);
np      = np-1;
cn      = approx.cofn;
cc      = approx.cofc;
bnds    = approx.bnds;
Ga1     = zeros(size(Ga,1),1);

I       = find(Ga<abar(iz));
J       = find(Ga>=abar(iz));

if ~isempty(I)
tmp     = 2*(log(Ga(I))-bnds(iz,1))/(bnds(iz,2)-bnds(iz,1))-1;
Ga1(I)= exp(cheb(tmp,np)*cn(:,iz));
end

if ~isempty(J)
tmp     = 2*(log(Ga(J))-bnds(iz,2))/(bnds(iz,3)-bnds(iz,2))-1;
Ga1(J)= exp(cheb(tmp,np)*cc(:,iz));
end