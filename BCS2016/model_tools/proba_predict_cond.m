function proba_predict_cond(file)
load(strcat('results/',file,'_simul'));
%%
th      = 0.1275;
T       = length(sim.y);
Dsc     = diff(sim.dum)==1; % look when a crisis starts in t+1
dum     = sim.dum(1:T-1);   % state in t (1-> crisis, 0 otherwise)
pr      = sim.pr1(1:T-1);   % P(crisis starts in t+1)
k       = sim.k(1:T-1);
y       = sim.y(1:T-1);
ksy     = sim.k(1:T-1)./sim.y(1:T-1);
R       = sim.R(1:T-1);
z       = sim.z(1:end-1);

Xm      = comp_type12(pr,dum,Dsc,th);
T0      = 3;
T       = length(pr);
pr      = max(pr(T0:T),1e-12);
x       = log(pr./(1-pr));

Dsc     = Dsc(T0:T);
dum     = dum(T0:T);
lk      = log(k(T0:T));
ly      = log(y(T0:T));
dk      = log(k(T0:T))-log(k(T0-1:T-1));
dy      = log(y(T0:T))-log(y(T0-1:T-1));
ksy     = log(ksy(T0:T));
R       = R(T0:T);
T       = length(x);
Dsc     = Dsc(dum==0);
x       = x(dum==0);
LY      = ly(dum==0);
LK      = lk(dum==0);
LZ      = z(dum==0);
LDK     = dk(dum==0);
LDY     = dy(dum==0);
LKSY    = ksy(dum==0);
LR      = R(dum==0);
cst     = ones(length(x),1);

%%
ROLS    = [NaN;NaN;Xm.*[100;100;1]];
res     = ols(x,[cst LZ]);
F       = (res.nobs-res.nvar)*res.rsqr/((res.nvar-1)*res.rsqr);
pF      = 1-fdis_cdf(F,res.nobs-res.nvar,res.nvar-1);
X       = comp_type12cond(1./(1+exp(-res.yhat)),Dsc,th);
ROLS    = [ROLS [res.rsqr;pF;X.*[100;100;1]]];

res     = ols(x,[cst LK]);
F       = (res.nobs-res.nvar)*res.rsqr/((res.nvar-1)*res.rsqr);
pF      = 1-fdis_cdf(F,res.nobs-res.nvar,res.nvar-1);
X       = comp_type12cond(1./(1+exp(-res.yhat)),Dsc,th);
ROLS    = [ROLS [res.rsqr;pF;X.*[100;100;1]]];

res     = ols(x,[cst LK LZ]);
F       = (res.nobs-res.nvar)*res.rsqr/((res.nvar-1)*res.rsqr);
pF      = 1-fdis_cdf(F,res.nobs-res.nvar,res.nvar-1);
X       = comp_type12cond(1./(1+exp(-res.yhat)),Dsc,th);
ROLS    = [ROLS [res.rsqr;pF;X.*[100;100;1]]];

res     = ols(x,[cst LKSY]);
F       = (res.nobs-res.nvar)*res.rsqr/((res.nvar-1)*res.rsqr);
pF      = 1-fdis_cdf(F,res.nobs-res.nvar,res.nvar-1);
X       = comp_type12cond(1./(1+exp(-res.yhat)),Dsc,th);
ROLS    = [ROLS [res.rsqr;pF;X.*[100;100;1]]];

XX      = [cst LKSY];
res     = logit(Dsc,XX);
pL      = 1-chis_prb(res.lratio,res.nvar-1);
X       = comp_type12cond(1./(1+exp(-XX*res.beta)),Dsc,th);
ROLS    = [ROLS [res.r2mf;pF;X.*[100;100;1]]];

%%
disp(' &Model &~& $k$ & $TFP$ & $(k,TFP)$ & $k/y$  &~&  $k/y$   \\')
str={'$R^2$','F-Test (\%)','Type I (\%)','Type II (\%)','Warnings'};
[r,c]   = size(ROLS);
for i=1:r-1
    fprintf(strcat('%20s ',repmat(' & %4.2f',1,c),' \\\\ \n'),str{i},ROLS(i,:));
end
fprintf(strcat('%20s ',repmat(' & %d',1,c),' \\\\ \n'),str{r},ROLS(r,:));

