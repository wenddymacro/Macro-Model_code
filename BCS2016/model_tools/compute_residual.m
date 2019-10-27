function compute_residual(stst,data,approx,param)

% Computes the phi function
% alpha   = param(1);
% gamma   = param(2);
% delta   = param(3);
nu      = param(4);
vth     = param(5);
% Rbar    = param(6);
% lambda  = param(7);
% theta   = param(8);
% pmin    = param(9);
psi     = param(10);
beta    = param(11);
sigma   = param(12);
amin    = exp(approx.bnds(1,1));
amax    = exp(approx.bnds(1,3));

abar    = data.abar;
Pz      = data.Pz;
Gz      = data.Gz;
nz      = length(Gz);
na      = 100;
at      = linspace(stst{1}.a1,stst{nz}.a1,na)';
dat.Ga  = at;
dat.Gz  = Gz;
eqt     = compute_equilibrium(dat,param);
income  = eqt.income;
h       = eqt.h;

res     = zeros(na,nz);
for i=1:nz;
    a1      = comp_ap(at,abar,approx,i);
    cons    = income(:,i)-psi*a1;
    uc      = cons-vth*h(:,i).^(1+nu)/(1+nu);
    a2      = comp_a2(a1,abar,approx);
    dat1.Ga = a1;
    dat1.Gz = Gz;
    eqtp1   = compute_equilibrium(dat1,param);
    r1      = eqtp1.r;
    income1 = eqtp1.income;
    h1      = eqtp1.h;
    expect  = (Pz(i,:)*((income1-psi*a2-vth*h1.^(1+nu)/(1+nu)).^(-sigma).*r1)')';
    res(:,i)= (uc-(beta*expect).^(-1/sigma))./cons;
end
disp('log10(E|res|) log10(E(res^2)) log10(max(|res|)')
E1  = log10(mean(abs(res)))';
E2  = log10(mean(res.^2))';
E3  = log10(max(abs(res)))';
for i=1:nz
    fprintf('%8.4f\t%8.4f\t%8.4f\n',E1(i),E2(i),E3(i))
end