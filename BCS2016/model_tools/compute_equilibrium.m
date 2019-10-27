function output=compute_equilibrium(data,param)

% Computes the phi function
alpha   = param(1);
gamma   = param(2);
delta   = param(3);
nu      = param(4);
vth     = param(5);
Rbar    = param(6);
lambda  = param(7);
theta   = param(8);
pmin    = param(9);


tol     = 1e-6;
Ga      = data.Ga;
Gz      = data.Gz;
na      = length(Ga);
nz      = length(Gz);

rho     = zeros(na,nz);
phi     = zeros(na,nz);
r       = zeros(na,nz);
R       = zeros(na,nz);
y       = zeros(na,nz);
k       = zeros(na,nz);
h       = zeros(na,nz);
income  = zeros(1,nz);
for j=1:nz
    
    z       = Gz(j);
    abar    = ((1-alpha)/vth)^(1/nu)*(alpha/(Rbar+delta-1))^((nu+alpha)/(nu*(1-alpha)))*z^((1+nu)/(nu*(1-alpha)));    
    htmp    = ((1-alpha)*z/vth).^(1/(alpha+nu))*Ga.^(alpha/(alpha+nu));
    Rtmp    = alpha*z*Ga.^(alpha-1).*htmp.^(1-alpha)+1-delta;

    % Normal times (a_t <= abar_t)
    J       = find(Ga<=abar);
    if ~isempty(J)
        nj      = length(J);
        k(J,j)  = Ga(J);
        h(J,j)  = ((1-alpha)*z/vth).^(1/(alpha+nu))*k(J,j).^(alpha/(alpha+nu));
        y(J,j)  = z*k(J,j).^alpha.*h(J,j).^(1-alpha);
        crit    = 1;
        pt0     = zeros(nj,1);
        while crit>tol
            p0  = (pmin+exp(pt0))./(1+exp(pt0));
            dp0 = (1-pmin)*exp(pt0)./((1+exp(pt0)).^2);
            f0  = Rtmp(J).*p0.*(1-p0.^lambda)+gamma*(1-theta)*(p0.^lambda)-gamma;
            df0 = (Rtmp(J).*(1-(lambda+1)*(p0.^lambda))+gamma*(1-theta)*lambda*(p0.^(lambda-1))).*dp0;
            pt1 = pt0-f0./df0;
            crit= max(abs(pt1-pt0));
            pt0 = pt1;
        end
        p           = (pmin+exp(pt0))./(1+exp(pt0));
        R(J,j)      = Rtmp(J);
        r(J,j)      = (lambda/(1+lambda))*Rtmp(J).*(1-p.^(lambda+1))./(1-p.^lambda);
        income(J,j) = y(J,j)+(1-delta)*k(J,j);
        rho(J,j)    = p.*Rtmp(J);
        phi(J,j)    = (rho(J,j)-gamma)/(gamma*theta);
    end;
    
    % Crisis times (a_t > abar_t)
    I       = find(Ga>abar);
    if ~isempty(I)
        crit    = 1;
        ka      = min(Ga)*ones(length(I),1);
        kb      = Ga(I);
        Rcst    = alpha*((1-alpha)/vth)^((1-alpha)/(alpha+nu))*z^((1+nu)/(alpha+nu));
        while crit>tol;
            k0  = (ka+kb)/2;
            Ra  = Rcst*ka.^(-nu*(1-alpha)/(alpha+nu))+1-delta;
            Rb  = Rcst*kb.^(-nu*(1-alpha)/(alpha+nu))+1-delta;
            R0  = Rcst*k0.^(-nu*(1-alpha)/(alpha+nu))+1-delta;
            fa  = Ra-gamma*(Ga(I)./(Ga(I)-ka)).^(1/lambda);
            fb  = Rb-gamma*(Ga(I)./(Ga(I)-kb)).^(1/lambda);
            f0  = R0-gamma*(Ga(I)./(Ga(I)-k0)).^(1/lambda);
            Ib  = sign(fa).*sign(f0)<0;
            Ia  = sign(fb).*sign(f0)<0;
            kb(Ib)  = k0(Ib);
            ka(Ia)  = k0(Ia);            
            crit= max(abs(f0));
        end
        k(I,j)      = k0;
        h(I,j)      = ((1-alpha)*z/vth).^(1/(alpha+nu))*k0.^(alpha/(alpha+nu));
        y(I,j)      = z*k0.^alpha.*h(I,j).^(1-alpha)+(gamma+delta-1)*(Ga(I)-k0);
        R(I,j)      = alpha*z*k0.^(alpha-1).*h(I,j).^(1-alpha)+(1-delta);
        r(I,j)      = (lambda+(gamma./R(I,j)).^(lambda+1)).*R(I,j)/(lambda+1);
        income(I,j) = y(I,j)+(1-delta)*Ga(I);
        rho(I,j)    = gamma;
        phi(I,j)    = 0;
    end
end

output.k        = k;
output.h        = h;
output.R        = R;
output.y        = y;
output.income   = income;
output.r        = r;
output.rho      = rho;
output.phi      = phi;