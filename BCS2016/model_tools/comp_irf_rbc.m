function comp_irf_rbc(file,varargin)
load(strcat('results/',file,'_approx'),'approx','stst','data','param')
if nargin>1
    filemc  = strcat('results/',varargin{1});
    simmc   = 0;
else
    simmc   = 1;
end
alpha   = param(1);
delta   = param(3);
nu      = param(4);
vth     = param(5);
psi     = param(10);
csth    = ((1-alpha)/vth)^(1/(alpha+nu));

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
    save('results/mc_irf','Z1p','Z2p');
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
zp1     = zeros(nsim,T);
hp1     = zeros(nsim,T);
cp1     = zeros(nsim,T);
yp1     = zeros(nsim,T);
kp1     = zeros(nsim,T);
xp1     = zeros(nsim,T);
zp2     = zeros(nsim,T);
hp2     = zeros(nsim,T);
cp2     = zeros(nsim,T);
yp2     = zeros(nsim,T);
kp2     = zeros(nsim,T);
xp2     = zeros(nsim,T);
for i=1:nsim
    kp1(i,1) = stst{i0}.k;
    kp2(i,1) = stst{i0}.k;
end
for t=1:T
    fprintf('period %d\n',t)
    zp1(:,t) = data.Gz(Z1p(:,t));
    hp1(:,t) = csth*zp1(:,t).^(1/(alpha+nu)).*kp1(:,t).^(alpha/(alpha+nu));
    yp1(:,t) = zp1(:,t).*kp1(:,t).^alpha.*hp1(:,t).^(1-alpha);
    
    zp2(:,t) = data.Gz(Z2p(:,t));
    hp2(:,t) = csth*zp2(:,t).^(1/(alpha+nu)).*kp2(:,t).^(alpha/(alpha+nu));
    yp2(:,t) = zp2(:,t).*kp2(:,t).^alpha.*hp2(:,t).^(1-alpha);
    for i=1:nz
        I           = find(Z1p(:,t)==i);
        kp1(I,t+1)  = comp_kp(kp1(I,t),approx,i);
        I           = find(Z2p(:,t)==i);
        kp2(I,t+1)  = comp_kp(kp2(I,t),approx,i);
    end
    xp1(:,t) = psi*kp1(:,t+1)-(1-delta)*kp1(:,t);
    cp1(:,t) = yp1(:,t)-xp1(:,t);
    xp2(:,t) = psi*kp2(:,t+1)-(1-delta)*kp2(:,t);
    cp2(:,t) = yp2(:,t)-xp2(:,t);
end
%%
yp  = 100*(log(yp1)-log(yp2));
hp  = 100*(log(hp1)-log(hp2));
cp  = 100*(log(cp1)-log(cp2));
xp  = 100*(log(xp1)-log(xp2));

my      = mean(yp)';
mc      = mean(cp)';
mx      = mean(xp)';
mh      = mean(hp)';
myrbc   = [my mh mc mx];
sy      = std(yp)';
sc      = std(cp)';
sx      = std(xp)';
sh      = std(hp)';
syrbc   = [sy sh sc sx];
save('results/irf_rbc','myrbc','syrbc');