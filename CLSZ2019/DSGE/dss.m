function d = dss(DELTA, G, OMEGAK, PSI, ETA, BETA, EPS, OMEGAP, PISS, ALPHA, THETA, MS,MP, XIS,XIP, K, OMEGAM, AS,AP, LS, LP, TAU,GSS,SIGMAM,PSIY)

% solve for the steady state

% solve for the steady state

hes = 1;
hep = 1;
I_K = 1-(1-DELTA)/G;% Eq. 4
rk = G/BETA-(1-DELTA); % Eq. 9; Eq. 8 gives uc = uk

r = G*PISS/BETA; % Eq. 7
x = EPS/(EPS-1); % Eq. 15
rs = (r-1)/(1-TAU)+1; % Eq. 31
rp = r; % Eq. 34
pi = PISS; % steady state inflation 
tau = TAU; % steady state RRR

% solve for bankruptcy cutoffs omegas omegap 
% using Eq. 16-28.
options = optimset('MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-16);
FF =fsolve(@(x)solveff(AS,AP,K,OMEGAM,LS,LP,MS,MP,XIS,XIP,ALPHA,THETA,G,PISS,BETA,TAU,OMEGAM*(1+exp(x(1))),OMEGAM*(1+exp(x(2)))),[log(1.02),log(1.01)],options);

omegas = OMEGAM*(1+exp(FF(1)));
omegap = OMEGAM*(1+exp(FF(2)));

FS = 1/(K-1)*OMEGAM^K*omegas^(-K+1); % fraction of revenue income going to the SOE borrower f(omegas)
GS = LS*omegas+(1-LS)*(1-(1-MS)*K/(K-1))*OMEGAM^K*omegas^(-K+1)+(1-MS)*(1-LS)*K/(K-1)*OMEGAM; % fraction of revenue income going to the SOE lender g(omegas)
G1S = LS+(1-LS)*(1-MS*K)*OMEGAM^K*omegas^(-K); % first order derivative of g(omegas)
F1S = -OMEGAM^K*omegas^(-K);% first order derivative of f(omegas)
BS_YS = GS/rs; % Eq. 26
NS_YS = -G1S/F1S*FS/rs*PISS*G; % Eq. 27
KS_YS = (1-ALPHA)*(NS_YS/PISS/G+BS_YS)/rk; % Eq. 22

FP = 1/(K-1)*OMEGAM^K*omegap^(-K+1); % fraction of revenue income going to the POE borrower f(omegap)
GP = LP*omegap+(1-LP)*(1-(1-MP)*K/(K-1))*OMEGAM^K*omegap^(-K+1)+(1-MP)*(1-LP)*K/(K-1)*OMEGAM;% fraction of revenue income going to the SOE lender g(omegap)
G1P = LP+(1-LP)*(1-MP*K)*OMEGAM^K*omegap^(-K);% first order derivative of g(omegap)
F1P = -OMEGAM^K*omegap^(-K);% first order derivative of g(omegap)
BP_YP = GP/rp; % Eq. 26
NP_YP = -G1P/F1P*FP/rp*PISS*G;% Eq. 27
KP_YP = (1-ALPHA)*(NP_YP/PISS/G+BP_YP)/rk; % Eq. 22

HSP1 =(NS_YS/PISS/G+BS_YS)/(NP_YP/PISS/G+BP_YP);
YSP = ((PSIY/(1-PSIY))^(SIGMAM/(SIGMAM-1))*AS/AP*(KS_YS/KP_YP)^(1-ALPHA)*HSP1^(ALPHA*THETA))^(1/(ALPHA-ALPHA*THETA+1/(SIGMAM-1)));%YSPP = YS/YP*PS/PP

%WSP = ((AS/AP)*(MU/(1-MU))^((1-THETA)*ALPHA/XIL)*(NS_YS/PISS/G+BS_YS)/(NP_YP/PISS/G+BP_YP))^(1/(ALPHA+(1-THETA)*ALPHA/XIL));
%YSP = WSP^((1-ALPHA+ALPHA*THETA)/XIL+1-ALPHA)/(MU/(1-MU))^((1-ALPHA+ALPHA*THETA)/XIL)*AS/AP;

YP_Y = 1/(1+YSP)/x;% YP_Y = YP*PP/(Y*P) = YP*PP/(YS*PS+YP*PP)*Pw/P
YS_Y = YSP/(1+YSP)/x;% YP_Y = YP*PP/(Y*P) = YP*PP/(YS*PS+YP*PP)*Pw/P
ps = (YS_Y/PSIY^SIGMAM/x^(-SIGMAM))^(1/(1-SIGMAM));
pp = (YP_Y/(1-PSIY)^SIGMAM/x^(-SIGMAM))^(1/(1-SIGMAM));

K_Y = (KS_YS*YS_Y+KP_YP*YP_Y)*G; % Eq. 46 
GDP_Y = 1-YS_Y*MS*K/(K-1)*(OMEGAM - OMEGAM^K*omegas^(-K+1))-YP_Y*MP*K/(K-1)*(OMEGAM - OMEGAM^K*omegap^(-K+1)); % Eq. 37
G_Y = GSS*GDP_Y; % Eq. 39
C_Y = GDP_Y-I_K*K_Y-G_Y; % Eq. 38

WAGE_Y = ALPHA*THETA*((NP_YP/PISS/G+BP_YP)*YP_Y+(NS_YS/PISS/G+BS_YS)*YS_Y); % Eq. 20
h = (WAGE_Y/C_Y/PSI)^(1/(ETA+1)); % Eq. 5,6
HS_HP = (NS_YS/PISS/G+BS_YS)/(NP_YP/PISS/G+BP_YP)*(YSP);
hhp = h/(1+(HS_HP));
hhs = h/(1+1/(HS_HP)); % Eq. 20 ???

ys = ((AS)*KS_YS^(1-ALPHA)*(hes^(1-THETA)*hhs^THETA)^ALPHA*ps^(1-ALPHA))^(1/ALPHA); % Eq. 16
yp = ((AP)*KP_YP^(1-ALPHA)*(hep^(1-THETA)*hhp^THETA)^ALPHA*pp^(1-ALPHA))^(1/ALPHA);% Eq. 16


y = ys/PSIY^SIGMAM/(ps*x)^(-SIGMAM);% Eq. 18
%y = yp/(1-PSIY)^SIGMAM/(pp*x)^(-SIGMAM);% Eq. 18

k = K_Y*y;
i = I_K*k;
c = C_Y*y;
g = G_Y*y;
uc = 1/c; % Eq. 5
w =PSI*h^(ETA)*c;
uk = uc;


WE_YP = ALPHA*(1-THETA)*(NP_YP/PISS/G+BP_YP);
WE_YS = ALPHA*(1-THETA)*(NS_YS/PISS/G+BS_YS);

wes =  WE_YS*ys*ps;
wep =  WE_YP*yp*pp;

ks = KS_YS*ys*ps;
ns = NS_YS*ys*ps;
bs = BS_YS*ys*ps;
as = ys*ps/(ns/PISS/G+bs);

kp = KP_YP*yp*pp;
np = NP_YP*yp*pp;
bp = BP_YP*yp*pp;
ap = yp*pp/(np/PISS/G+bp);

d = (bs)/(1-tau) + (bp);

end



