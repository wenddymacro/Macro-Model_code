function plot_typicpath(file)
load(strcat('results/',file,'_approx'));   % loads the file containing approximation data
load(strcat('results/',file,'_simul'));    % loads the file containing simulated data
pause_  = 1;                    % To pause the graphs
nz      = length(Gz);
ss      = stst{(nz+1)/2};
%%
hor     = length(sim.y);

tb      = 25;
te      = 10;
T       = tb+te+1;

tc      = find(diff(sim.dum)==1)+1;         % Find initial date of crisis
tc      = tc(and((tc<=hor-te),(tc>=tb)));   % Defines sample
Nc      = length(tc);                       % Number of crises
ma      = zeros(Nc,T);
mla     = zeros(Nc,T);
me      = zeros(Nc,T);
mab     = zeros(Nc,T);
mz      = zeros(Nc,T);
mk      = zeros(Nc,T);
mh      = zeros(Nc,T);
mc      = zeros(Nc,T);
mx      = zeros(Nc,T);
my      = zeros(Nc,T);
mR      = zeros(Nc,T);
mr      = zeros(Nc,T);
mphi    = zeros(Nc,T);
mda     = zeros(Nc,T);
msp     = zeros(Nc,T);
mrho    = zeros(Nc,T);
msiz    = zeros(Nc,T);
mpr1    = zeros(Nc,T);
maphi   = zeros(Nc,T);
mks     = zeros(Nc,T);
mdum    = zeros(Nc,T);
miz     = zeros(Nc,T);
mlevl   = zeros(Nc,T);
mlevb   = zeros(Nc,T);
mleva   = zeros(Nc,T);
mky     = zeros(Nc,T);
msr     = zeros(Nc,T);
ee      = filter([1 -0.9],1,log(sim.z));

for i=1:Nc;
    smpl        = tc(i)-tb:tc(i)+te;
    mz(i,:)     = log(sim.z(smpl));
    miz(i,:)    = sim.iz(smpl);
    ma(i,:)     = sim.a(smpl);
    mla(i,:)    = log(sim.a(smpl));
    me(i,:)     = ee(smpl);
    mk(i,:)     = sim.k(smpl);
    my(i,:)     = sim.y(smpl);
    mh(i,:)     = sim.h(smpl);
    mx(i,:)     = sim.x(smpl);
    mc(i,:)     = sim.c(smpl);
    mR(i,:)     = sim.R(smpl);
    mr(i,:)     = sim.r(smpl);
    mab(i,:)    = sim.abar(smpl);
    msp(i,:)    = sim.spread(smpl);
    mphi(i,:)   = sim.phi(smpl);
    mrho(i,:)   = sim.rho(smpl);
    msiz(i,:)   = sim.bsize(smpl);
    mpr1(i,:)   = sim.pr1(smpl);
    maphi(i,:)  = sim.aggphi(smpl);
    mda(i,:)    = 100*(sim.a(smpl+1)-sim.a(smpl))./sim.y(smpl);
    mks(i,:)    = sim.k(smpl)./sim.a(smpl);
    mky(i,:)    = sim.k(smpl)./sim.y(smpl);
    mdum(i,:)   = sim.dum(smpl);    
    mlevl(i,:)  = sim.levl(smpl);
    mlevb(i,:)  = sim.levb(smpl);
    mleva(i,:)  = sim.leva(smpl);
    msr(i,:)    = sim.sr(smpl);
end

%%
smpl    = 1:T;
Psit    = 1.^(smpl-1)';
zt      = zeros(T,1);
for t=smpl
   zt(t)= mean(mdum(:,t))>0.5;
end
pct     = [100/6 500/6];
step    = 5:8:T;
ms      = 4;
mt2     = 's';
op      ='median';
xt      =(-tb:te)';
%%
y       = ma(:,smpl);
subplot(221);plotexpband(xt,[feval(op,y)'.*Psit prctile(y,pct)'.*repmat(Psit,1,2)],[ss.abar.*Psit],zt,'Assets & Absorption Capacity');hold on
plot(xt,ss.a1.*Psit,'linewidth',2,'color','k','linestyle','--') 
plot(xt,feval(op,mab(:,smpl))'.*Psit,'linewidth',1,'color','k','linestyle','-')
set(gca,'xlim',[-tb te],'ylim',[2.5 4.75])
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)

y       = mz(:,smpl);
subplot(222);plotexpband(xt,[feval(op,y)'+log(Psit) prctile(y,pct)'+repmat(log(Psit),1,2)],[0],zt,'TFP Level (Log.)');hold on
set(gca,'xlim',[-tb te],'ylim',[-0.1 0.1])
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y1=median(me(:,smpl))';
plot(xt,y1,'linewidth',1,'color','k','linestyle','-')
plot(xt(step),y1(step)',mt2,'markersize',ms,'markeredgecolor','k','markerfacecolor','k')

% if pause_;pause;end
print('-depsc2','figures/typicpath1');
close
%%
fs2=10;
y=mR(:,smpl);
subplot(3,4,1);plotexpband(xt,[feval(op,y)' prctile(y,pct)'],100*([ss.R]-1),zt,'Corporate Loan Rate');ylabel('percents','fontname','times','fontsize',10)
set(gca,'xlim',[-tb te],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y=mr(:,smpl);
subplot(3,4,2);plotexpband(xt,[feval(op,y)' prctile(y,pct)'],100*([ss.r]-1),zt,'Return on Deposits');ylabel('percents','fontname','times','fontsize',10)
set(gca,'xlim',[-tb te],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y=mrho(:,smpl);
subplot(3,4,3);plotexpband(xt,[feval(op,y)' prctile(y,pct)'],100*([ss.rho]-1),zt,'Interbank Rate');ylabel('percents','fontname','times','fontsize',10)
set(gca,'xlim',[-tb te],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y=mphi(:,smpl);
subplot(3,4,4);plotexpband(xt,[feval(op,y)' prctile(y,pct)'],[ss.phi],zt,{'Market Funding Ratio','(Borrowing Bank)'});
set(gca,'xlim',[-tb te],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y=msiz(:,smpl);
subplot(3,4,5);plotexpband(xt,[feval(op,y)'.*Psit prctile(y,pct)'.*repmat(Psit,1,2)],[ss.bsize*Psit],zt,'Size of Banking Sector');
set(gca,'xlim',[-tb te],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y=mks(:,smpl);
subplot(3,4,6);plotexpband(xt,[feval(op,y)' prctile(y,pct)'],[ss.k]./[ss.a1],zt,'Credit/Assets');
set(gca,'xlim',[-tb te],'ylim',[0.875 1.025],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y=mky(:,smpl);
subplot(3,4,7);plotexpband(xt,[feval(op,y)' prctile(y,pct)'],[ss.k]./[ss.y],zt,'Credit/Output');
set(gca,'xlim',[-tb te],'ylim',[1.9 2.4],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y=mpr1(:,smpl);
subplot(3,4,8);plotexpband(xt,[feval(op,y)' prctile(y,pct)'],[ss.pr1],zt,'1-step ahead Proba.');
set(gca,'xlim',[-tb te],'ylim',[-0.05 1],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y=mh(:,smpl);
subplot(3,4,9);plotexpband(xt,[feval(op,y)' prctile(y,pct)'],[ss.h],zt,'Hours Worked');
set(gca,'xlim',[-tb te],'ylim',[0.9 1.3],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y       = my(:,smpl);
subplot(3,4,10);plotexpband(xt,[feval(op,y)'.*Psit prctile(y,pct)'.*repmat(Psit,1,2)],[ss.y.*Psit],zt,'Output');
set(gca,'xlim',[-tb te],'ylim',[1.2 2],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y=mc(:,smpl);
subplot(3,4,11);plotexpband(xt,[feval(op,y)'.*Psit prctile(y,pct)'.*repmat(Psit,1,2)],[ss.c.*Psit],zt,'Consumption');
set(gca,'xlim',[-tb te],'ylim',[0.9 1.5],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
y=mx(:,smpl);
subplot(3,4,12);plotexpband(xt,[feval(op,y)'.*Psit prctile(y,pct)'.*repmat(Psit,1,2)],[ss.x.*Psit],zt,'Investment');
set(gca,'xlim',[-tb te],'ylim',[0.25 0.55],'fontsize',fs2)
h=line([0 0],get(gca,'ylim'));set(h,'color','k','linestyle','-.','linewidth',1)
% if pause_;pause;end
print('-depsc2','figures/typicpath2');
close
