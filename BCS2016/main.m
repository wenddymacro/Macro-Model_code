% This is the main file to produce most of the results for the paper
% "Booms and Banking Crises" by F. Boissay, F. Collard and F. Smets
%
% It is advisable to run the file main_rbc first.
%
clear all;clc;close all
format short
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
if ~exist('figures','dir')
    mkdir('figures')
end
if ~exist('results/rbc_approx.mat','file')
    error('Please run main_rbc.m before running this file.')
end
file    = 'benchmark';  % Name of results files
solve_  = 1;            % 1-> solves the model, 0-> loads computed solution
simul_  = 1;            % 1-> simulates the model, 0-> loads existing simulation
%%
% decentralized economy: structural parameters

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
amax    = 8;            % maximum value of savings

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
param   = [alpha gamma delta nu vth Rbar lambda theta pmin psi beta sigma];
disp('Deterministic Steady State:')
fprintf('R: %4.2f%%\tr: %4.2f%%\tgamma: %4.2f%%\tR-r: %4.2f%%\n',100*(Rss-1),100*(rss-1),100*(gamma-1),spread)
fprintf('p: %g\tpmin: %g\tk/y: %4.2f\tK: %4.2f\t\tKbar: %4.2f\n',pss,pmin,kss/yss,kss,kbar)

%% Stochastic Process and Associated Threshold
mth         = 0;    % 1 -> Rouwenhorst, 0 -> Tauchen-Hussey
[Gz,Pz]     = build_markov_chain(rz,sez,nz,mth);
abar        = ((1-alpha)/vth)^(1/nu)*(alpha/(Rbar+delta-1))^((nu+alpha)/(nu*(1-alpha)))*Gz.^((1+nu)/(nu*(1-alpha)));
data.Pz     = Pz;
data.Gz     = Gz;
data.abar   = abar;
data.kss    = kss;
data.tol    = 1e-6;
disp(' ')
disp('    Z         abar')
disp([Gz(:) abar(:)])

%% Grid for savings
zx      = rcheb(na);
bnds    = zeros(nz,3);
Ga      = zeros(2*na,nz);
inc     = zeros(2*na,nz);
h       = zeros(2*na,nz);
for i=1:nz;
    bnds(i,:)       = [log(amin) log(abar(i)) log(amax)];
    Ga(1:na,i)      = exp((1+zx)*(bnds(i,2)-bnds(i,1))/2+bnds(i,1));
    Ga(na+1:2*na,i) = exp((1+zx)*(bnds(i,3)-bnds(i,2))/2+bnds(i,2));
    dat1.Ga         = Ga(:,i);
    dat1.Gz         = Gz(i);
    eq_t            = compute_equilibrium(dat1,param);
    inc(:,i)        = eq_t.income;
    h(:,i)          = eq_t.h;
end
data.Ga = Ga;
data.inc= inc;
data.h  = h;
%% Initial Values
Gag     = repmat(linspace(amin,amax,200)',1,nz);
Ga1     = 0.05+0.9*Ga;
approx  = smoothdecrule(Ga,Ga1,abar,amin,amax,np);
%% Main Loop
switch solve_
    case 1
        dat1.Gz = Gz;
        crit    = 1;
        iter    = 1;
        while crit>tol
            Ga1     = comp_a1(Ga,abar,approx);              % Compute a(t+1)
            hh      = plot(Gag,Gag,'k--',Gag,comp_a1(Gag,abar,approx));
            set(hh,'linewidth',1);pause(0.01)
            
            a10     = zeros(2*na,nz);
            for j=1:nz;
                dat1.Ga = Ga1(:,j);
                Ga2     = comp_a2(Ga1(:,j),abar,approx);    % Compute all possible a(t+2)
                equil   = compute_equilibrium(dat1,param);  % Compute Equilibrium in t+1
                r1      = equil.r;                          % Obtain r(t+1)
                inc1    = equil.income;                     % Obtain income(t+1)
                h1      = equil.h;                          % Obtain h(t+1)
                expect  = (Pz(j,:)*((inc1-psi*Ga2-vth*h1.^(1+nu)/(1+nu)).^(-sigma).*r1)')';
                a10(:,j)= (inc(:,j)-vth*h(:,j).^(1+nu)/(1+nu)-(beta*expect).^(-1/sigma))/psi;
            end
            tmp         = smoothdecrule(Ga,a10,abar,amin,amax,np);            
            crit        = max(max(abs([tmp.cofn;tmp.cofc]-[approx.cofn;approx.cofc])));            
            approx.cofn = learn*tmp.cofn+(1-learn)*approx.cofn;
            approx.cofc = learn*tmp.cofc+(1-learn)*approx.cofc;            
            fprintf('iteration #%5d\tcriterion = %g\n',iter,crit)
            iter        = iter+1;
        end
        %% Compute Euler residuals
        %     compute_residual(stst,data,approx,param)
        
        %% Compute riskfree rate
        apprf   = compute_riskfree_rate(approx,data,param,np);
        
        %% Compute Steady State reached in each TFP level
        stst    = compute_steady_state(approx,apprf,data,param);
        
        %% Some bookkeeping
        save(strcat('results/',file,'_approx'),'approx','apprf','Ga','Gz','Pz','stst','param','data');
        
    otherwise
        load(strcat('results/',file,'_approx'))
end

%% Simulations & moments
switch simul_
    case 1
        data.T      = 500000;               % Length of data set
        data.burn   = 1000;                 % Burning period
        data.initz  = (size(Pz,2)+1)/2;     % Initial state (starts from average)
        sim         = simulate_dynamics(approx,apprf,data,param,stst,1);%,'results/markovchain');
        save(strcat('results/',file,'_simul'),'sim')
    otherwise
        load(strcat('results/',file,'_simul'))
end
pr          = mean(diff(sim.dum)>0);
fprintf('Frequncy of crises: %4.2f%%\n',100*pr)
disp('    rho       R         r         R-rf')
disp([mean(sim.rho(sim.dum==0)) mean([sim.R sim.r sim.R-sim.rf])])

%% Impulse response functions to a one-standard deviation shock
comp_irf(file,'mc_irf')

%% Plotting Typical path
plot_typicpath(file);

%% Statistics
table   = 1;    % 1-> generates latex table, 0-> just reports results
gr_     = 1;    % 1-> generates figures (model vs data), 0-> no graphs
comp_stats(file,table,gr_);

%% Crisis Prediction Analysis
proba_predict_cond(file)