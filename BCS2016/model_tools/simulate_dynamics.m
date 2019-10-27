function out = simulate_dynamics(approx,apprf,data,param,stst,varargin)
%
% Computes simulated data based on the solution decision rules
% s0 is the initial level of savings
% ib is the initial state for the Markov chain
% data contains:
%   data.Ga : grid for a_t
%   data.Pa : Transition matrix for a_t
%   data.T  : # of periods
%   data.burn: # of periods for burn out stage
%
T       = data.T+data.burn;
burn    = data.burn;
smpl    = burn+1:T;
nz      = size(data.Pz,2);
if nargin>6
    load(varargin{2})
    if length(iz)<T
        error('Markov Chain is too short')
    end
else
    iz      = simulate_markov(1:nz,data.Pz,T,data.initz);
end
out.iz  = iz(smpl);
gamma   = param(2);
lambda  = param(7);
a       = zeros(T+1,1);
k       = zeros(T,1);
x       = zeros(T,1);
c       = zeros(T,1);
y       = zeros(T,1);
h       = zeros(T,1);
R       = zeros(T,1);
Rt      = zeros(T,1);
Rf      = zeros(T,1);
r       = zeros(T,1);
r       = zeros(T,1);
p       = zeros(T,1);
phi     = zeros(T,1);
z       = zeros(T,1);
dum     = zeros(T,1);
rho     = zeros(T,1);
sp      = zeros(T,1);
ab      = zeros(T,1);
pr1     = zeros(T,1);
aggphi  = zeros(T,1);
bsz     = zeros(T,1);
ncore   = zeros(T,1);
core    = zeros(T,1);
cash    = zeros(T,1);
dse     = zeros(T,1);
sr      = zeros(T,1);
if nargin>5
    a(1)    = stst{data.initz}.a1;
    for t=1:T
        if rem(t,25000)==0;fprintf('%5d periods computed\n',t);end
        if nargin>5
            if varargin{1}==0
                output  = compute_dynamics(a(t),iz(t),data,approx,apprf,param);
            else
                output  = compute_dynamics(a(t),iz(t),data,approx,apprf,param,1);
                pr1(t)  = output.pr1;
            end
        end
        a(t+1)  = output.a1;
        z(t)    = output.z;
        ab(t)   = output.abar;
        
        k(t)    = output.k;
        h(t)    = output.h;
        c(t)    = output.c;
        x(t)    = output.x;
        y(t)    = output.y;
        
        p(t)    = output.p;
        R(t)    = output.R;
        Rt(t)   = output.Rtmp;
        Rf(t)   = output.rf;
        r(t)    = output.r;
        sp(t)   = output.spread;
        rho(t)  = output.rho;
        phi(t)  = output.phi;
        dum(t)  = output.dum;
        bsz(t)  = output.bsize;
        ncore(t)= output.ncore;
        core(t) = output.core;
        cash(t) = output.cashhold;
        aggphi(t)= output.aggphi;
        dse(t)  = gamma/max(Rf(t)-gamma,0);
        sr(t)   = output.sr;
    end
    
    out.z       = z(smpl);
    out.a       = a(smpl);
    out.abar    = ab(smpl);
    
    out.k       = k(smpl);
    out.h       = h(smpl);
    out.c       = c(smpl);
    out.x       = x(smpl);
    out.y       = y(smpl);
    
    out.R       = 100*(R(smpl)-1);
    out.Rt      = 100*(Rt(smpl)-1);
    out.r       = 100*(r(smpl)-1);
    out.rf      = 100*(Rf(smpl)-1);
    out.phi     = phi(smpl);
    out.rho     = 100*(rho(smpl)-1);
    out.dum     = dum(smpl);
    out.spread  = 100*sp(smpl);
    out.bsize   = bsz(smpl);
    out.ncore   = ncore(smpl);
    out.core    = core(smpl);
    out.cash    = cash(smpl);
    out.aggphi  = aggphi(smpl);
    out.dse     = dse(smpl-1);
    out.levl    = dse(smpl-1);
    out.levb    = (1-dum(smpl)).*(dse(smpl-1)+phi(smpl).*(1+dse(smpl-1)))+dum(smpl).*dse(smpl-1);
    out.leva    = (1-dum(smpl)).*(p(smpl).^lambda.*out.levl+(1-p(smpl).^lambda).*out.levb)+dum(smpl).*out.levl;
    out.sr      = sr(smpl);
    if nargin>5
        if varargin{1}~=0
            out.pr1 = pr1(smpl);
        end
    end
    
else
    a(1)    = stst{data.initz}.a1;
    for t=1:T
        if rem(t,25000)==0;fprintf('%5d periods computed\n',t);end
        output  = compute_dynamics(a(t),iz(t),data,approx,apprf,param,1);
        a(t+1)  = output.a1;
        z(t)    = output.z;
        ab(t)   = output.abar;
        
        k(t)    = output.k;
        h(t)    = output.h;
        c(t)    = output.c;
        x(t)    = output.x;
        y(t)    = output.y;
        
        p(t)    = output.p;
        R(t)    = output.R;
        Rt(t)   = output.Rtmp;
        Rf(t)   = output.rf;
        r(t)    = output.r;
        sp(t)   = output.spread;
        rho(t)  = output.rho;
        phi(t)  = output.phi;
        dum(t)  = output.dum;
        bsz(t)  = output.bsize;
        ncore(t)= output.ncore;
        core(t) = output.core;
        cash(t) = output.cashhold;
        aggphi(t)= output.aggphi;
        dse(t)  = gamma/max(Rf(t)-gamma,0);
        pr1(t)  = output.pr1;
        sr(t)   = output.sr;
    end
    smpl        = burn+1:T;
    out.z       = z(smpl);
    out.a       = a(smpl);
    out.abar    = ab(smpl);
    
    out.k       = k(smpl);
    out.h       = h(smpl);
    out.c       = c(smpl);
    out.x       = x(smpl);
    out.y       = y(smpl);
    
    out.R       = 100*(R(smpl)-1);
    out.Rt      = 100*(Rt(smpl)-1);
    out.r       = 100*(r(smpl)-1);
    out.rf      = 100*(Rf(smpl)-1);
    out.phi     = phi(smpl);
    out.rho     = 100*(rho(smpl)-1);
    out.dum     = dum(smpl);
    out.spread  = 100*sp(smpl);
    out.bsize   = bsz(smpl);
    out.ncore   = ncore(smpl);
    out.core    = core(smpl);
    out.cash    = cash(smpl);
    out.aggphi  = aggphi(smpl);
    out.dse     = dse(smpl-1);
    out.levl    = dse(smpl-1);
    out.levb    = (1-dum(smpl)).*(dse(smpl-1)+phi(smpl).*(1+dse(smpl-1)))+dum(smpl).*dse(smpl-1);
    out.leva    = (1-dum(smpl)).*(p(smpl).^lambda.*out.levl+(1-p(smpl).^lambda).*out.levb)+dum(smpl).*out.levl;
    out.sr      = sr(smpl);
    
    out.pr1     = pr1(smpl);
end