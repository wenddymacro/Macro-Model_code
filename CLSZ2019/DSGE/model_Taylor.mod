%% dynare file to solve for simple rule under stochastic simulation 
%% with both fiscal shock and productivity shock
%% respond to GDP gap
var k i uc c h w r pi uk rk x y wes wep 
ys ks hes hhs ns bs as rs omegas
yp kp hep hhp np bp ap rp omegap
z_s tau gdp ytfp gtfp rzs rzp utility totutility g bcs bcp q levs levp omegas1 omegap1 ps pp a_p a_s diffa m m2 d m2g bps yps ltau;


varexo  eps_s eps_g eps_ap eps_as eps_da;

parameters DELTA G OMEGAK PSI ETA BETA EPS OMEGAP PISS ALPHA THETA  XIS XIP K OMEGAM AP AS LS LP TAU RP RY GSS MS MP TY TP RHOS RHOG RR TT PSIY SIGMAM GDPSS SIGS SIGG SIGAP RHOAP SIGAS RHOAS SIGDA RHODA PSIM DSS;

parameter=calibration1;
DELTA=parameter(1);
G=parameter(2);
OMEGAK=parameter(3);
PSI=parameter(4);
PSIY=parameter(5);
ETA=parameter(6);
BETA=parameter(7);
EPS=parameter(8);
OMEGAP=parameter(9);
PISS=parameter(10);
ALPHA=parameter(11);
THETA=parameter(12);
XIS=parameter(13);
XIP=parameter(14);
K=parameter(15);
OMEGAM=parameter(16);
AS=parameter(17);
AP=parameter(18);
LS=parameter(19);
LP=parameter(20);
TAU=parameter(21);
GSS=parameter(22);
SIGMAM=parameter(23);
MS=parameter(24);
MP=parameter(25);
SIGS=parameter(26);
SIGG=parameter(27);
SIGAP=parameter(28);
SIGAS=parameter(29);
SIGDA=parameter(30);

GDPSS = gdpss(DELTA, G, OMEGAK, PSI, ETA, BETA, EPS, OMEGAP, PISS, ALPHA, THETA, MS,MP, XIS,XIP, K, OMEGAM, AS,AP, LS, LP, TAU,GSS,SIGMAM,PSIY);
CSS   = css(DELTA, G, OMEGAK, PSI, ETA, BETA, EPS, OMEGAP, PISS, ALPHA, THETA, MS,MP, XIS,XIP, K, OMEGAM, AS,AP, LS, LP, TAU,GSS,SIGMAM,PSIY);
DSS   = dss(DELTA, G, OMEGAK, PSI, ETA, BETA, EPS, OMEGAP, PISS, ALPHA, THETA, MS,MP, XIS,XIP, K, OMEGAM, AS,AP, LS, LP, TAU,GSS,SIGMAM,PSIY);

PSIM  = (1 - BETA/PISS/G)/(CSS)*(DSS/0.95*0.05);%currency stands for 5% of M2

RHOS = 0.95;
RHOG = 0.95;
RHOAP = 0.95;
RHOAS = 0.95;
RHODA = 0.95;

RR = 0;
RP = -0.397/(1-0.391);
RY = 0.183/(1-0.391);
TP = 0; 
TY = 0; 
TT = 0;


model ;

%% Household sector Equation 
exp(h) = exp(hhs)+exp(hhp);
exp(k) = (1-DELTA)*exp(k(-1))/G+(1-OMEGAK/2*(exp(i-i(-1))*G-G)^2)*exp(i);
uc =-c;
exp(uc+w) = PSI*exp((ETA)*h); 
exp(uc) = BETA*exp(r+uc(+1)-pi(+1)-ln(G));
exp(uc) = exp(uk)*(1-OMEGAK/2*(exp(i-i(-1))*G-G)^2-OMEGAK*(exp(i-i(-1))*G-G)*exp(i-i(-1))*G)+BETA*exp(uk(+1))/G*OMEGAK*(exp(i(+1)-i)*G-G)*(exp(i(+1)-i)*G)^2;
exp(uk) = BETA*(exp(uk(+1))/G*(1-DELTA)+exp(uc(+1)+rk(+1))/G);
q = uk-uc;
exp(uc) = BETA*exp(uc(+1)-pi(+1)-ln(G)) + PSIM*exp(-m);

%% Retail 
exp(-x) = (EPS-1)/EPS + OMEGAP/EPS*exp(-y)*((exp(pi)/PISS-1)*exp(pi+c)/PISS-BETA*exp(uc(+1)-uc)*(exp(pi(+1))/PISS-1)*exp(pi(+1)+c(+1))/PISS);
%% Wholesale
exp(y*(SIGMAM-1)/SIGMAM) = (PSIY*exp(ys*(SIGMAM-1)/SIGMAM)+(1-PSIY)*exp(yp*(SIGMAM-1)/SIGMAM));
ys = SIGMAM*ln(PSIY)-SIGMAM*(ps+x)+y;
yp = SIGMAM*ln(1-PSIY)-SIGMAM*(pp+x)+y;

%% Firms Equation 18,22-24,27-29
ys = log(AS)+z_s-diffa/2+a_s+(1-ALPHA)*ks+ALPHA*((1-THETA)*hes+THETA*hhs);
exp(w+hhs) = ALPHA*THETA*(exp(ns(-1)-pi)/G+exp(bs));
exp(wes+hes)= ALPHA*(1-THETA)*(exp(ns(-1)-pi)/G+exp(bs));
exp(rk+ks)= (1-ALPHA)*(exp(ns(-1)-pi)/G+exp(bs));
exp(as) = exp(ys+ps)/(exp(ns(-1)-pi)/G+exp(bs));
exp(as)*(exp(ns(-1)-pi)/G+exp(bs))*(LS*exp(omegas)+(1-LS)*(1-(1-MS)*K/(K-1))*OMEGAM^K*exp(omegas)^(-K+1)+(1-MS)*(1-LS)*K/(K-1)*OMEGAM) = exp(rs+bs);
exp(ns(-1)-pi)/G/(exp(ns(-1)-pi)/G+exp(bs)) = (LS+(1-LS)*(1-MS*K)*OMEGAM^K*exp(omegas)^(-K))/(OMEGAM^K*exp(omegas)^(-K))*exp(as)*(1/(K-1)*OMEGAM^K*exp(omegas)^(-K+1))/exp(rs);
exp(ns) = XIS*exp(as)*(exp(ns(-1)-pi)/G+exp(bs))*(1/(K-1)*OMEGAM^K*exp(omegas)^(-K+1))+exp(wes);


yp =log(AP)+z_s+diffa/2+a_p+(1-ALPHA)*kp+ALPHA*((1-THETA)*hep+THETA*hhp);
exp(w+hhp) = ALPHA*THETA*(exp(np(-1)-pi)/G+exp(bp));
exp(wep+hep)= ALPHA*(1-THETA)*(exp(np(-1)-pi)/G+exp(bp));
exp(rk+kp)= (1-ALPHA)*(exp(np(-1)-pi)/G+exp(bp));
exp(ap) = exp(yp+pp)/(exp(np(-1)-pi)/G+exp(bp));
exp(ap)*(exp(np(-1)-pi)/G+exp(bp))*(LP*exp(omegap)+(1-LP)*(1-(1-MP)*K/(K-1))*OMEGAM^K*exp(omegap)^(-K+1)+(1-MP)*(1-LP)*K/(K-1)*OMEGAM) = exp(rp+bp);
exp(np(-1)-pi)/G/(exp(np(-1)-pi)/G+exp(bp)) = (LP+(1-LP)*(1-MP*K)*OMEGAM^K*exp(omegap)^(-K))/(OMEGAM^K*exp(omegap)^(-K))*exp(ap)*(1/(K-1)*OMEGAM^K*exp(omegap)^(-K+1))/exp(rp);
exp(np) = XIP*exp(ap)*(exp(np(-1)-pi)/G+exp(bp))*(1/(K-1)*OMEGAM^K*exp(omegap)^(-K+1))+exp(wep);

%% Financial intermediaries Equation 32 35
(exp(rs)-1)*(1-(tau)) = exp(r)-1;
rp = r;

%% monetary policy 36 37
tau - (1-TT)*TAU =TT*tau(-1)+TY*(gdp-gdp(-1)) +TP*(pi-ln(PISS));

exp(d) = exp(bs)/(1-tau) + exp(bp);
exp(m2) =  exp(m) + exp(bs)/(1-tau);

m2g = m2 - m2(-1) + pi + ln(G);
m2g - (1-RR)*ln(G*PISS) = RR*m2g(-1)+(RY*(gdp-gdp(-1)) + RP* (pi -ln(PISS)));


%% Market clearing 38-47 2
exp(y) = exp(c) + OMEGAP/2*(exp(pi)/PISS-1)^2*exp(c) + exp(i) + exp(g)+ exp(as)*(exp(ns(-1)-pi)/G+exp(bs))*MS*K/(K-1)*(OMEGAM - OMEGAM^K*exp(omegas)^(-K+1)) + exp(ap)*(exp(np(-1)-pi)/G+exp(bp))*MP*K/(K-1)*(OMEGAM - OMEGAM^K*exp(omegap)^(-K+1));
exp(gdp) = exp(c) +  exp(i)+exp(g) ;
0 = hes;
0 = hep;
exp(k(-1))/G = exp(ks) + exp(kp);

%% Shock process 
z_s = RHOS*z_s(-1)+eps_s; %TFP shock
(g-gdp) - ln(GSS) = RHOG*(g(-1)-gdp(-1)-ln(GSS)) + eps_g;

a_p = RHOAP*a_p(-1)+eps_ap; %POE TFP shock
a_s = RHOAS*a_s(-1)+eps_as; %SOE TFP shock
diffa = RHODA*diffa(-1)+eps_da;% difference in POE and SOE

%% utilities and TFP
exp(ytfp) = (exp(y))/((exp(ks)+exp(kp))^(1-ALPHA)*(exp(h)^(THETA))^ALPHA);
exp(gtfp) = exp(gdp)/((exp(ks)+exp(kp))^(1-ALPHA)*(exp(h)^(THETA))^ALPHA);
utility = (c-PSI*exp(h)^(1+ETA)/(1+ETA));

totutility = utility + BETA*totutility(+1);
%% calculate SOEs' and POEs' loan rate with risk premium
exp(rzs)  = exp(omegas+as - bs)*(exp(ns(-1)-pi)/G+exp(bs));
exp(rzp)  = exp(omegap+ap - bp)*(exp(np(-1)-pi)/G+exp(bp));
%% calculate SOEs' and POEs' bankruptch ratio
exp(bcp) = 1-(OMEGAM/exp(omegap))^K;
exp(bcs) = 1-(OMEGAM/exp(omegas))^K;
exp(omegas) = (OMEGAM*(1+exp(omegas1)));
exp(omegap) = (OMEGAM*(1+exp(omegap1)));
levs = exp(bs - (ns(-1)-pi - log(G)));
levp = exp(bp - (np(-1)-pi - log(G)));

yps = exp(yp + pp - ys - ps);
bps = exp(bp - bs);
ltau = log(tau);
end;


steady_state_model ;
[k, i, uc, c, h ,w, r ,pi, uk, rk, x, y, wes,wep, ys, ks, hes, hhs, ns, bs, as ,rs, omegas, yp, kp, hep, hhp, np, bp, ap, rp, omegap, tau,gdp,ytfp,gtfp,rzs,rzp,g,bcs,bcp,ps,pp] = solvesslog(DELTA, G, OMEGAK, PSI, ETA, BETA, EPS, OMEGAP, PISS, ALPHA, THETA, MS,MP, XIS,XIP, K, OMEGAM, AS,AP, LS, LP, TAU,GSS,SIGMAM,PSIY);
m = (c-log((1 - BETA/PISS/G)/PSIM));
//utility = (c-PSI*exp(h)^(1+ETA)/(1+ETA) + PSIM*(m));
utility = (c-PSI*exp(h)^(1+ETA)/(1+ETA));

q = 0;
levs = exp(bs - (ns-pi-log(G)));
levp = exp(bp - (np-pi-log(G)));
totutility = utility/(1-BETA);
omegas1 = log(exp(omegas)/OMEGAM-1);
omegap1 = log(exp(omegap)/OMEGAM-1);

d = log(exp(bs)/(1-tau) + exp(bp));
m2 =  log(exp(m) + exp(bs)/(1-tau));
//m2 = log(exp(m) + tau*exp(bs)/(1-tau));
m2g = pi + log(G);
bps = exp(bp - bs);
yps = exp(yp + pp - ys - ps);
ltau = log(tau);

end;

shocks;
var eps_s;stderr SIGS; 
var eps_g;stderr SIGG;
var eps_ap;stderr SIGAP;
var eps_as;stderr SIGAS;
var eps_da;stderr SIGDA;
end;

resid;
steady;
check;
model_diagnostics;
stoch_simul(nograph,order=2,irf = 20);


