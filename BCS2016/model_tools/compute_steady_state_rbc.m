function as=compute_steady_state(approx,data,param)
% as=compute_steady_state(approx,apprf,data,param)
% 
% Input:    approx  : Structure for approximation coefficients of a(t+1)
%           apprf   : Structure for approximation coefficients of riskfree rate
%           data    : Structure that passes information for the decision
%                     rules (see main file)
%           param   : Parameters
%
% Output:   as      : Cell structure of nz elements, each containing the
%                     steady state associated with each level of the technology shock
%           Example : as{12}.y contains the steady state output level for the
%                     12th value of the technology shock.
%
%% Computes steady states
alpha   = param(1);
delta   = param(3);
nu      = param(4);
vth     = param(5);
psi     = param(10);
csth    = ((1-alpha)/vth)^(1/(alpha+nu));

lb      = 0.5;
a0      = data.kss;
itmax   = 2000;
tol     = data.tol;
Gz      = data.Gz;
nz      = size(Gz,1);
as      = cell(nz,1);
for i=1:nz
    crit    = 1;
    iter    = 1;
    while and(crit>tol,iter<=itmax)
        a1  = comp_kp(a0,approx,i);
        crit= abs(a1-a0);
        a0  = (1-lb)*a0+lb*a1;
        iter= iter+1;
    end
    as{i}.h = csth*Gz(i)^(1/(alpha+nu))*a0.^(alpha/(alpha+nu));
    as{i}.y = Gz(i)*a0.^alpha.*as{i}.h.^(1-alpha);
    as{i}.k = a0;
    as{i}.x = (psi+delta-1)*a0;
    as{i}.c = as{i}.y-as{i}.x;
end
