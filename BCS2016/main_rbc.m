clear all;clc;close all;
addpath('gen_tools')
addpath('model_tools')
%%
rng(1234567890,'twister')
learn   = 1;
tol     = 1e-6;
np      = 15;
na      = np+1;
nz      = 15;
if ~exist('results','dir')
    mkdir('results')
end
file    = 'rbc';        % Name of result file
solve_  = 1;            % set to 1 to solve the model, 0 just to load existing results.
%%
% decentralized economy, model with ex-ante deadweight loss
% structural parameters

% YEARLY PARAMETRIZATION
beta    = 1/1.03;       % Discount rate (adjusted for growth)
alpha   = 0.3;          % Capital elasticity
hss     = 1;            % h in det. steady state
nu      = 0.5;          % Inverse Frish labor elasticity
delta   = 0.1;          % Depreciation rate
sigma   = 4.5;          % Inverse intertemp. elasticity of subst.
Rbar    = 1.0262;       % Threshold Corporate loan rate
lambda  = 25;           % Distribution parameter
theta   = 0.093;        % Diversion parameter
rz      = 0.9;          % Persistence of TFP shock
sez     = 0.0177;       % Std. Dev. of TFP innovation
psi     = 1.012;        % Gross rate of growth
amin    = 0.5;          % minimum value of savings
amax    = 14;           % maximum value of savings

pmin    = ((2+(lambda-1)*theta-sqrt((2+(lambda-1)*theta)^2-4*(1-theta)))/(2*(1-theta)))^(1/lambda);
gamma   = Rbar*pmin*(1-pmin^lambda)/(1-(1-theta)*pmin^lambda);
vp      = (pmin:0.0000001:0.9999999)';
R1      = gamma*(1-(1-theta)*vp.^lambda)./(vp.*(1-vp.^lambda));
R2      = (1+lambda)*(1-vp.^lambda)./(beta*lambda*(1-vp.^(lambda+1)));
[~,ind] = min(abs(R1-R2));
pss     = vp(ind);
Rss     = R1(ind);
rss     = 1/beta;
spread  = 100*(Rss-rss);
kss     = (alpha*beta/(Rss+delta-1))^(1/(1-alpha));
yss     = kss^alpha;
ksy     = kss/yss;
kbar    = (alpha/(Rbar+delta-1))^(1/(1-alpha));
vth     = (1-alpha)*kss^alpha;
csth    = ((1-alpha)/vth)^(1/(alpha+nu));
param   = [alpha gamma delta nu vth Rbar lambda theta pmin psi beta sigma];

disp('Deterministic Steady State:')
fprintf('R: %4.2f%%\tr: %4.2f%%\tgamma: %4.2f%%\tR-r: %4.2f%%\n',100*(Rss-1),100*(rss-1),100*(gamma-1),spread)
fprintf('p: %g\tpmin: %g\tk/y: %4.2f\tK: %4.2f\t\tKbar: %4.2f\n',pss,pmin,kss/yss,kss,kbar)

%% Stochastic Process and Associated Threshold
mth         = 0;    % 1 -> Rouwenhorst, 0 -> Tauchen-Hussey
[Gz,Pz]     = build_markov_chain(rz,sez,nz,mth);

data.Pz     = Pz;
data.Gz     = Gz;
data.kss    = kss;
data.tol    = 1e-6;

%% Grid for savings
zx      = rcheb(na);
bnds    = zeros(nz,2);
Ga      = zeros(na,nz);
y       = zeros(na,nz);
h       = zeros(na,nz);
for i=1:nz;
    bnds(i,:)   = [log(amin) log(amax)];
    Ga(1:na,i)  = exp((1+zx)*(bnds(i,2)-bnds(i,1))/2+bnds(i,1));
    dat1.Ga     = Ga(:,i);
    dat1.Gz     = Gz(i);
    h(:,i)      = csth*Gz(i)^(1/(alpha+nu))*Ga(:,i).^(alpha/(alpha+nu));
    y(:,i)      = Gz(i)*Ga(:,i).^alpha.*h(:,i).^(1-alpha);
end
data.Ga = Ga;
data.Gz = Gz;
data.Pz = Pz;
%% Initial Values
Gag     = repmat(linspace(amin,amax,200)',1,nz);
Ga1     = 0.05+0.9*Ga;
approx  = smoothdecrule_rbc(Ga,Ga1,amin,amax,np);

%% Main Loop
if solve_
    dat1.Gz = Gz;
    crit    = 1;
    iter    = 1;
    while crit>tol
        Ga1     = comp_k1(Ga,approx);
        
        a10     = zeros(na,nz);
        hh      = plot(Gag,comp_k1(Gag,approx),Gag,Gag,'k--');
        pause(0.01)
        for j=1:nz;
            Ga2     = comp_k2(Ga1(:,j),approx);
            h1      = zeros(na,nz);
            y1      = zeros(na,nz);
            r1      = zeros(na,nz);
            c1      = zeros(na,nz);
            for i=1:nz
                h1(:,i)  = csth*Gz(i)^(1/(alpha+nu))*Ga1(:,j).^(alpha/(alpha+nu));
                y1(:,i)  = Gz(i)*Ga1(:,j).^alpha.*h1(:,i).^(1-alpha);
                r1(:,i)  = alpha*y1(:,i)./Ga1(:,j)+1-delta;
                c1(:,i)  = y1(:,i)+(1-delta)*Ga1(:,j)-psi*Ga2(:,i);
            end
            
            expect  = (Pz(j,:)*((c1-vth*h1.^(1+nu)/(1+nu)).^(-sigma).*r1)')';
            a10(:,j)= (y(:,j)+(1-delta)*Ga(:,j)-vth*h(:,j).^(1+nu)/(1+nu)-(beta*expect).^(-1/sigma))/psi;
        end
        approx1 = smoothdecrule_rbc(Ga,a10,amin,amax,np);
        crit    = max(max(abs(approx1.coef-approx.coef)));
        approx.coef   = learn*approx1.coef+(1-learn)*approx.coef;
        fprintf('iteration #%5d\tcriterion = %g\n',iter,crit)
        iter    = iter+1;
    end
    
    %% Compute Steady State
    stst    = compute_steady_state_rbc(approx,data,param);
    
    %% Some bookkeeping
    save(strcat('results/',file,'_approx'),'approx','Ga','Gz','Pz','stst','param','data');
    
else
    load(strcat('results/',file,'_approx'))
end
comp_irf_rbc(file)