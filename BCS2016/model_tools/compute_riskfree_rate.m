function apprf=compute_riskfree_rate(approx,data,param,np)

% Parameters
alpha   = param(1);
gamma   = param(2);
delta   = param(3);
nu      = param(4);
vth     = param(5);
Rbar    = param(6);
lambda  = param(7);
theta   = param(8);
pmin    = param(9);
psi     = param(10);
beta    = param(11);
sigma   = param(12);

% Data
Ga      = data.Ga;
Gz      = data.Gz;
Pz      = data.Pz;
abar    = data.abar;
[na,nz] = size(Ga);
% Computation of the riskfree rate
Ga1     = comp_a1(Ga,abar,approx);
muc     = (data.inc-psi*Ga1-vth*data.h.^(1+nu)/(1+nu)).^(-sigma);
rf      = zeros(na,nz);

for j=1:nz;
    Ga2     = comp_a2(Ga1(:,j),abar,approx);
    dat1.Ga = Ga1(:,j);
    dat1.Gz = Gz;
    equil   = compute_equilibrium(dat1,param);
    inc1    = equil.income;
    h1      = equil.h;
    muc1    = (inc1-psi*Ga2-vth*h1.^(1+nu)/(1+nu)).^(-sigma);
    expect  = (Pz(j,:)*muc1')';
    rf(:,j) = muc(:,j)./(beta*expect);
end

% Approximation
amin    = exp(approx.bnds(1,1));
amax    = exp(approx.bnds(1,3));
apprf   = smoothdecrule(Ga,rf,abar,amin,amax,np);