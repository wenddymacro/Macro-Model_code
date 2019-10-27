function comp_irf(file,varargin)
load(strcat('results/',file,'_approx'),'approx','stst','data','param')
if nargin>1
    filemc  = strcat('results/',varargin{1});
    simmc   = 0;
else
    simmc   = 1;
end
disp(' ')
disp(' Computation of IRFs')
disp('=====================')
disp(' ')

nz      = size(data.Pz,1);
i0      = (nz+1)/2;
if simmc
    nsim    = 100000;
    T       = 40;
    Z1p     = zeros(nsim,T);
    Z2p     = zeros(nsim,T);
    disp('    * Generating Markov Chains')
    disp('      ========================')
    start   = tic;
    for s   = 1:nsim
        [Z1p(s,:),Z2p(s,:)]	= simulate_markov2(1:nz,data.Pz,T,[i0+1;i0]);
    end
    timelap = toc(start);
    fprintf('It took %g sec. to generate %d Markov chains\n\n',timelap,nsim)
else
    load(filemc);
    [nsim,T]= size(Z1p);
end
%%
%
% shock
%
disp('    * Generating Dynamics')
disp('      ===================')
zp      = zeros(nsim,T);
ap1     = zeros(nsim,T+1);
ap2     = zeros(nsim,T+1);
hp      = zeros(nsim,T);
cp      = zeros(nsim,T);
yp      = zeros(nsim,T);
Rp      = zeros(nsim,T);
rp      = zeros(nsim,T);
kp      = zeros(nsim,T);
xp      = zeros(nsim,T);
spreadp = zeros(nsim,T);
phip    = zeros(nsim,T);
rhop    = zeros(nsim,T);
sizp    = zeros(nsim,T);
pr1p    = zeros(nsim,T);
aggphip = zeros(nsim,T);

ap1(:,1)= stst{i0}.a1;
ap2(:,1)= stst{i0}.a1;
for t=1:T
    fprintf('period %d\n',t)
    outputp1= compute_vec_dynamics(ap1(:,t),Z1p(:,t),data,approx,param,1);
    outputp2= compute_vec_dynamics(ap2(:,t),Z2p(:,t),data,approx,param,1);
    
    zp(:,t)     = 100*(log(outputp1.z)-log(outputp2.z));
    ap1(:,t+1)  = outputp1.a1;
    ap2(:,t+1)  = outputp2.a1;
    kp(:,t)     = outputp1.k./ap1(:,t);
    hp(:,t)     = 100*(log(outputp1.h)-log(outputp2.h));
    cp(:,t)     = 100*(log(outputp1.c)-log(outputp2.c));
    xp(:,t)     = 100*(log(outputp1.x)-log(outputp2.x));
    yp(:,t)     = 100*(log(outputp1.y)-log(outputp2.y));
    spreadp(:,t)= 100*outputp1.spread;
    rhop(:,t)   = 100*(outputp1.rho-1);
    Rp(:,t)     = 100*(outputp1.R-1);
    rp(:,t)     = 100*(outputp1.r-1);
    phip(:,t)   = outputp1.phi;
    sizp(:,t)   = outputp1.bsize;
    pr1p(:,t)   = outputp1.pr1;
    aggphip(:,t)= outputp1.aggphi;
end
%%
load('results/irf_rbc') % myrbc contains y,h,c,x
fs  = 10;
step= 5:5:40;

i   = 1;
my  = mean(yp)';
sy  = std(yp)';
subplot(241);h=plotirf([my myrbc(:,i)],sy);
title('Output','fontname','times','fontsize',fs)
ylabel('% deviation','fontname','times','fontsize',fs)
set(h(2),'color','k','linewidth',1,'linestyle','-')
hold on
h=plot(step,myrbc(step,i),'o');
set(h,'markeredgecolor','k','markerfacecolor','k','markersize',4)

i   = 3;
my  = mean(cp)';
sy  = std(cp)';
subplot(242);h=plotirf([my myrbc(:,i)],sy);
title('Consumption','fontname','times','fontsize',fs)
ylabel('% deviation','fontname','times','fontsize',fs)
set(h(2),'color','k','linewidth',1,'linestyle','-')
hold on
h=plot(step,myrbc(step,i),'o');
set(h,'markeredgecolor','k','markerfacecolor','k','markersize',4)

i   = 4;
my  = mean(xp)';
sy  = std(xp)';
subplot(243);h=plotirf([my myrbc(:,i)],sy);
title('Investment','fontname','times','fontsize',fs)
ylabel('% deviation','fontname','times','fontsize',fs)
set(h(2),'color','k','linewidth',1,'linestyle','-')
hold on
h=plot(step,myrbc(step,i),'o');
set(h,'markeredgecolor','k','markerfacecolor','k','markersize',4)

i   = 2;
my  = mean(hp)';
sy  = std(hp)';
subplot(244);h=plotirf([my myrbc(:,i)],sy);
title('Hours Worked','fontname','times','fontsize',fs)
ylabel('% deviation','fontname','times','fontsize',fs)
set(h(2),'color','k','linewidth',1,'linestyle','-')
hold on
h=plot(step,myrbc(step,i),'o');
set(h,'markeredgecolor','k','markerfacecolor','k','markersize',4)

% pause
print('-depsc2','figures/irfs');
close
