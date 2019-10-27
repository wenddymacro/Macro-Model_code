function output=compute_dynamics(a,iz,data,approx,apprf,param,varargin)
%
% Computes the equilibrium allocation for a pair (a_t,z_t)
% a is the level of savings
% iz is the state for the Markov chain
% data must contain:
%   data.abar   : Threshold values for a
%   data.Gz     : grid for z_t
%   data.Pz     : Transition matrix for z_t
%
% approx must contain:
%   approx.bnds : Bounds for the Chebychev transformation
%   approx.cofn : coefficients in normal time
%   approx.cofc : coefficients during the crisis


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
% beta    = param(11);
sigma   = param(12);

tol     = 1e-6;
abar    = data.abar;
Gz      = data.Gz;
Pz      = data.Pz;
z       = Gz(iz);
a1      = comp_ap(a,abar,approx,iz);
rf      = comp_ap(a,abar,apprf,iz);
htmp    = ((1-alpha)*z/vth).^(1/(alpha+nu)).*a.^(alpha/(alpha+nu));
Rtmp    = alpha*z.*a.^(alpha-1).*htmp.^(1-alpha)+1-delta;

if a<abar(iz)
    % Normal times
    k        = a;
    h        = ((1-alpha)*z/vth).^(1/(alpha+nu)).*k.^(alpha/(alpha+nu));
    y        = z.*k.^alpha.*h.^(1-alpha);
    crit     = 1;
    pt0      = 0;
    while crit>tol
        p0  = (pmin+exp(pt0))./(1+exp(pt0));
        dp0 = (1-pmin)*exp(pt0)./((1+exp(pt0)).^2);
        f0  = Rtmp.*p0.*(1-p0.^lambda)+gamma*(1-theta)*(p0.^lambda)-gamma;
        df0 = (Rtmp.*(1-(lambda+1)*(p0.^lambda))+gamma*(1-theta)*lambda*(p0.^(lambda-1))).*dp0;
        pt1 = pt0-f0./df0;
        crit= max(abs(pt1-pt0));
        pt0 = pt1;
    end
    p       = (pmin+exp(pt0))./(1+exp(pt0));
    R       = Rtmp;
    r       = (lambda/(1+lambda))*Rtmp.*(1-p.^(lambda+1))./(1-p.^lambda);
    income  = y+(1-delta)*k;
    conso   = income-psi*a1;
    x       = income-(1-delta)*a-conso;
    dum      = 0;
    % leverage and interbank rate
    rho      = p.*R;
    phi      = (rho-gamma)/(gamma*theta);
    aggphi  = (1-p.^lambda).*phi;
    % size of the banking sector as a whole
    bsize    = a.*(1+p.^lambda);
    % share of core and non-core assets/liabilities and cash holding
    ncore    = (1-p.^lambda).*phi.*a; % same as a.*-p.^lambda
    core     = bsize-ncore;
    cashhold = 0;
else
    % crisis times
    k0      = a;
    crit    = 1;
    Rcst    = alpha*((1-alpha)/vth)^((1-alpha)/(alpha+nu))*z^((1+nu)/(alpha+nu));
    while crit>tol;
        R0  = Rcst*k0.^(-nu*(1-alpha)/(alpha+nu))+1-delta;
        dR0 = -nu*(1-alpha)*(R0+delta-1)./((alpha+nu)*k0);
        f0  = (a-k0).*R0.^lambda-(gamma^lambda)*a;
        df0 = -R0.^lambda+lambda*(a-k0).*(R0.^(lambda-1)).*dR0;
        k1  = k0-f0./df0;
        crit= max(abs(k1-k0));
        k0  = k1;
    end
    k       = k0;
    h       = ((1-alpha)*z/vth).^(1/(alpha+nu))*k0.^(alpha/(alpha+nu));
    y       = z*k0.^alpha.*h.^(1-alpha)+(gamma+delta-1)*(a-k0);
    R       = alpha*z*k0.^(alpha-1).*h.^(1-alpha)+(1-delta);
    r       = (lambda+(gamma./R).^(lambda+1)).*R/(lambda+1);
    income  = y+(1-delta)*a;
    conso   = income-psi*a1;
    x       = income-(1-delta)*a-conso;
    dum      = 1;
    % leverage and interbank rate
    phi      = 0;
    rho      = gamma;        % here by convention
    p        = gamma./R;  % here by convention
    aggphi  = (1-p.^lambda).*phi;
    % size of the banking sector
    bsize    = a;
    % share of core and non-core assets/liabilities and cash holding
    ncore  = 0;
    core   = a;
    cashhold = p.^lambda;
end

if nargin>6
    % Probability that next crisis breaks out in t+1
    % Find the threshold abar for at+1 (given st+1)
    ztmp        = (vth/(1-alpha))^((1-alpha)/(1+nu))*((Rbar+delta-1)/alpha)^((nu+alpha)/(1+nu));
    ztmp1       = ztmp*a1.^(nu*(1-alpha)/(1+nu));
    output.pr1  = Pz(iz,:)*(Gz<ztmp1);
%     output.pr2  = 0;
%     a2          = zeros(1,nz);
%     for i2=1:nz
%         a2(i2)  = comp_ap(a1,abar,approx,i2);
%     end
%     % Probabilities that next crisis breaks out in t+2, t+3, t+4
%     iz1     = find(Gz>ztmp1);
%     n1      = length(iz1);
%     for j1  = 1:n1
%         % for all the at+1>atmp1 and given st+1 find st+2(st+1,at+1)
%         % Probability of a crisis
%         ztmp2       = ztmp*a2(iz1(j1))^(nu*(1-alpha)/(1+nu));
%         output.pr2  = output.pr2+ Pz(iz,iz1(j1))*Pz(iz1(j1),:)*(Gz<=ztmp2);
%     end
end
output.z        = z;
output.y        = y;
output.k        = k;
output.x        = x;
output.h        = h;
output.c        = conso;
output.a1       = a1;
output.abar     = abar(iz);
output.dum      = dum;
output.sr       = 1-conso./(conso+x);    % Saving rate
output.income   = income;
output.uc       = (conso-vth*h^(1+nu)/(1+nu))^(-sigma);

output.R        = R;
output.r        = r;
output.rf       = rf;
output.rho      = rho;
output.spread   = R-r;
output.Rtmp     = Rtmp;

output.p        = p;
output.phi      = phi;
output.aggphi   = aggphi;
output.bsize    = bsize;
output.ncore  = ncore;
output.core   = core;
output.cashhold = cashhold;
