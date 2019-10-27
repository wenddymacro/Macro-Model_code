function output=compute_vec_dynamics(a,iz,data,approx,param,varargin)
%
% Computes the equilibrium allocation for a pair (s_t,a_t)
% s is the level of savings
% ia is the state for the Markov chain
%
% data must contain:
%   data.Gs     : grid for s_t
%   data.Ga     : grid for a_t
%   data.Pa     : Transition matrix for a_t
%   data.sbar   : Threshold value for s
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
a       = a(:);
na      = length(a);


abar    = data.abar;
Gz      = data.Gz;
Pz      = data.Pz;
nz      = length(Gz);
z       = Gz(iz);
z       = z(:);
ab      = abar(iz);
ab      = ab(:);
a1      = zeros(na,1);
for i=1:nz;
    I       = find(iz==i);
    if ~isempty(I)
        a1(I)   = comp_ap(a(I),abar,approx,i);
    end
end

htmp    = ((1-alpha)*z/vth).^(1/(alpha+nu)).*a.^(alpha/(alpha+nu));
Rtmp    = alpha*z.*a.^(alpha-1).*htmp.^(1-alpha)+1-delta;

I       = find((a-ab)<=0);
ni      = length(I);
J       = find((a-ab)>0);
nj      = length(J);

k       = zeros(na,1);
x       = zeros(na,1);
h       = zeros(na,1);
y       = zeros(na,1);
p       = zeros(na,1);
r       = zeros(na,1);
rho     = zeros(na,1);
income  = zeros(na,1);
phi     = zeros(na,1);
bsize   = zeros(na,1);
aggphi  = zeros(na,1);
R       = zeros(na,1);
shcore  = zeros(na,1);
shncore = zeros(na,1);
cashhold= zeros(na,1);
dum     = zeros(na,1);
if ni>0
    % Normal times
    k(I)        = a(I);
    h(I)        = ((1-alpha)*z(I)/vth).^(1/(alpha+nu)).*k(I).^(alpha/(alpha+nu));
    y(I)        = z(I).*k(I).^alpha.*h(I).^(1-alpha);
    crit        = 1;
    pt0         = zeros(ni,1);
    while crit>tol
        p0  = (pmin+exp(pt0))./(1+exp(pt0));
        dp0 = (1-pmin)*exp(pt0)./((1+exp(pt0)).^2);
        f0  = Rtmp(I).*p0.*(1-p0.^lambda)+gamma*(1-theta)*(p0.^lambda)-gamma;
        df0 = (Rtmp(I).*(1-(lambda+1)*(p0.^lambda))+gamma*(1-theta)*lambda*(p0.^(lambda-1))).*dp0;
        pt1 = pt0-f0./df0;
        crit= max(abs(pt1-pt0));
        pt0 = pt1;
    end
    p(I)     = (pmin+exp(pt0))./(1+exp(pt0));
    R(I)     = Rtmp(I);
    r(I)     = (lambda/(1+lambda))*Rtmp(I).*(1-p(I).^(lambda+1))./(1-p(I).^lambda);
    income(I)= y(I)+(1-delta)*k(I);
    dum(I)   = 0;
    % leverage and interbank rate
    rho(I)   = p(I).*R(I);
    phi(I)   = (rho(I)-gamma)/(theta*gamma);
    aggphi(I)= (1-p(I).^lambda).*phi(I);
    % size of the banking sector as a whole
    bsize(I) = a(I).*(1+p(I).^lambda);
    % share of core and non-core assets/liabilities and cash holding
    shncore(I)  = (1-p(I).^lambda).*phi(I)./bsize(I);
    shcore(I)   = 1-shncore(I);
    cashhold(I) = 0;
    
end

if nj>0
    % crisis times
    Rcst    = alpha*((1-alpha)/vth)^((1-alpha)/(alpha+nu)).*z(J).^((1+nu)/(alpha+nu));
    k0      = a(J);
    crit    = 1;
    while crit>tol;
        R0  = Rcst.*k0.^(-nu*(1-alpha)/(alpha+nu))+1-delta;
        dR0 = -nu*(1-alpha)*(R0+delta-1)./((alpha+nu)*k0);
        f0  = (a(J)-k0).*R0.^lambda-(gamma^lambda)*a(J);
        df0 = -R0.^lambda+lambda*(a(J)-k0).*(R0.^(lambda-1)).*dR0;
        k1  = k0-f0./df0;
        crit= max(abs(k1-k0));
        k0  = k1;
    end
    k(J)        = k0;
    h(J)        = ((1-alpha)*z(J)/vth).^(1/(alpha+nu)).*k0.^(alpha/(alpha+nu));
    y(J)        = z(J).*k0.^alpha.*h(J).^(1-alpha)+(gamma+delta-1)*(a(J)-k0);
    R(J)        = alpha*z(J).*k0.^(alpha-1).*h(J).^(1-alpha)+1-delta;
    r(J)        = (lambda+(gamma./R(J)).^(lambda+1)).*R(J)/(lambda+1);
    income(J)   = y(J)+(1-delta)*a(J);
    dum(J)      = 1;
    % leverage and interbank rate
    phi(J)      = 0;
    rho(J)      = gamma;        % here by convention
    p(J)        = gamma./R(J);  % here by convention
    aggphi(J)   = (1-p(J).^lambda).*phi(J);
    % size of the banking sector
    bsize(J)    = a(J);
    % share of core and non-core assets/liabilities and cash holding
    shncore(J)  = 0;
    shcore(J)   = 1;
    cashhold(J) = p(J).^lambda;
end

if nargin>5
    % Probability that next crisis breaks out in t+1
    % Find the threshold abar for at+1 (given st+1)
    ztmp        = (vth/(1-alpha))^((1-alpha)/(1+nu))*((Rbar+delta-1)/alpha)^((nu+alpha)/(1+nu));
    ztmp1       = ztmp*a1.^(nu*(1-alpha)/(1+nu));
    d1          = zeros(na,nz);
    for i=1:nz
        d1(:,i) = (ztmp1>=Gz(i));
    end
    tmp=sum(Pz(iz,:).*d1,2);
    output.pr1  = tmp;
%     output.pr2  = zeros(na,1);
%     for i=1:na
%         % Probabilities that next crisis breaks out in t+2, t+3, t+4
%         iz1         = find(Gz>ztmp1(i));
%         n1          = length(iz1);
%         for j1= 1:n1
%             % for all the at+1>atmp1 and given st+1 find st+2(st+1,at+1)
%             % Probability of a crisis
%             ztmp2           = ztmp*a2(i,iz1(j1))^(nu*(1-alpha)/(1+nu));
%             output.pr2(i)   = output.pr2(i)+ Pz(iz(i),iz1(j1))*Pz(iz1(j1),:)*(Gz<=ztmp2);
%         end
%     end    
end
output.a1       = a1;
output.k        = k;
output.h        = h;
output.R        = R;
output.y        = y;
output.c        = income-psi*a1;
output.x        = income-(1-delta)*a-output.c;
output.sr       = output.x./(output.c+output.x);
output.phi      = phi;
output.aggphi   = aggphi;
output.r        = r;
output.abar     = ab;
output.spread   = R-r;
output.z        = z;
output.Rtmp     = Rtmp;
output.dum      = dum;
output.p        = p;
output.income   = income;
output.rho      = rho;
output.bsize    = bsize;
output.shncore  = shncore;
output.shcore   = shcore;
output.cashhold = cashhold;