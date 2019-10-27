close all

var Cs Hs Ns D Rs Ws q Cb Hb Nb B Rb Wb k Cf lamb CRR lamf Y A
dY dB SPREAD dCs dCb dCf dHs dHb dq;
varexo eps_A;

%predetermined_variables B; 

parameters betas jj etaa betab betaf alphaa rhoA sigmaA kSS CRRSS
CsSS HsSS NsSS DSS RsSS WsSS qSS CbSS HbSS NbSS BSS RbSS WbSS CfSS lambSS lamfSS YSS ASS
gammaa SIG SIG2 SIG3;

betas =0.99;
betab =0.98;
jj     =0.1;
etaa   =2;
betaf =0.965;
alphaa =0.64;
rhoA  =0.9;
sigmaA=0.01;
kSS   =0.8992;
CRRSS =0.001;

paramsaux = [betas betab jj etaa betaf alphaa rhoA sigmaA kSS CRRSS];
start = [1.04093];

options = optimset('Display','iter','algorithm','Levenberg-Marquardt','FunValCheck','off','MaxFunEvals',500000,'MaxIter',500000);
[solution,fval] = fsolve(@(x) calcss(x,paramsaux),start,options);

NbSS = solution(1);

gammaa = 1-CRRSS;
SIG = (1-betab-(1-betab/betaf*(1-(1-betaf*1/betas)*gammaa))*kSS*((1-(1-betaf*1/betas)*gammaa)*1/betaf)^(-1));
SIG2 = (1/betas*gammaa+1-gammaa-(1-(1-betaf*1/betas)*gammaa)*1/betaf)*kSS*((1-(1-betaf*1/betas)*gammaa)*1/betaf)^(-1);
SIG3 = -gammaa*kSS*jj*(1-1/betas)*(1-alphaa)/alphaa*SIG^(-1)*((1-(1-betaf*1/betas)*gammaa)*1/betaf)^(-1);
%NbSS = (1+kSS*jj*SIG^(-1)*(((1-(1-betaf*1/betas)*gammaa)*1/betaf)^(-1)-1))^(1/etaa);
NsSS = (SIG3*NbSS^(-etaa)+1)^(1/-etaa);
CfSS = -SIG2*jj*NsSS^(alphaa)*NbSS^(1-alphaa-etaa)*(1-alphaa)*SIG^(-1);
HbSS = (alphaa/(1-alphaa)*1/(1-betas)*SIG*NsSS^(-etaa)*NbSS^(etaa)+1)^(-1);
qSS = jj/(1-HbSS)*1/(1-betas)*(alphaa*NsSS^(alphaa-etaa)*NbSS^(1-alphaa));
RbSS = (1-(1-betaf*1/betas)*gammaa)*1/betaf;
CbSS = (1-alphaa)*NsSS^(alphaa)*NbSS^(-alphaa-etaa+1);
CsSS = alphaa*NsSS^(alphaa-etaa)*NbSS^(1-alphaa);
lambSS = 1/CbSS*(1-RbSS*betab);
lamfSS = 1/CfSS*(1-betaf*1/betas);
HsSS = 1-HbSS;
BSS = 1/RbSS*kSS*qSS*HbSS;
DSS = gammaa*BSS;
RsSS = 1/betas;
WsSS = NsSS^(etaa-1)*CsSS;
WbSS = NbSS^(etaa-1)*CbSS;
YSS = NsSS^(alphaa)*NbSS^(1-alphaa);
ASS = 1;

model;
%Savers
1/Cs = betas*(1/Cs(+1)*Rs);%1
q/Cs = jj/Hs+betas*q(+1)/Cs(+1);%2
Ws = Ns^(etaa-1)*Cs;%3

%Borrowers
lamb*(B - 1/Rb(+1)*k*q(+1)*Hb)=0;%4
1/Cb = betab*1/Cb(+1)*Rb+lamb;%5
jj/Hb = 1/Cb*q-betab*q(+1)/Cb(+1)-lamb*1/Rb(+1)*k*q(+1);%6
Wb = Nb^(etaa-1)*Cb;%7

%Financial Intermediaries
lamf*(D - (1-CRR)*B)=0;%8
1/Cf = betaf*1/Cf(+1)*Rs+lamf;%9
1/Cf = betaf*1/Cf(+1)*Rb(+1)+(1-CRR)*lamf;%10

%Captial Requirement Ratio
CRR = CRRSS;%*(B/B(-1))^(0.1)*(Y/YSS)^(1.9)*(q/qSS)^(1.6);%11

%Loan-to-Value Ratio
k = kSS;%12

%Firms
Y = A*Ns^(alphaa)*Nb^(1-alphaa);%13
Ws = alphaa*Y/Ns;%14
Wb = (1-alphaa)*Y/Nb;%15

%Equilibrium
%Y = Cs + Cb + Cf;%16
Hs + Hb = 1;%16
Cf + Rs(-1)*D(-1)+B = D + Rb*B(-1);%17
Cb + Rb*B(-1) + q*(Hb-Hb(-1)) = B + Wb*Nb;%18
Cs + D + q*(Hs-Hs(-1)) = Rs(-1)*D(-1) + Ws*Ns;%19

%Exogenous Process
A = A(-1)^(rhoA)*exp(eps_A);%20


%Plot
dY = (Y-YSS)/YSS*100;
dB = (B-BSS)/BSS*100;
SPREAD = (Rb-Rs)*100;
dCs = (Cs-CsSS)/CsSS*100;
dCb = (Cb-CbSS)/CbSS*100;
dCf = (Cf-CfSS)/CfSS*100;
dHs = (Hs-HsSS)/HsSS*100;
dHb = (Hb-HbSS)/HbSS*100;
dq = (q-qSS)/qSS*100;


end;

steady_state_model;
Cs = CsSS;
Hs = HsSS;
Ns = NsSS;
D  = DSS;
Rs = RsSS;
Ws = WsSS;
q  = qSS;
Cb = CbSS;
Hb = HbSS;
Nb = NbSS;
B  = BSS;
Rb = RbSS;
Wb = WbSS;
k  = kSS;
Cf = CfSS;
lamb = lambSS;
CRR  = CRRSS;
lamf = lamfSS;
Y = YSS;
A = 1;
dY = 0;
dB = 0;
SPREAD = (RbSS-RsSS)*100;
dCs = 0;
dCb = 0;
dCf = 0;
dHs = 0;
dHb = 0;
dq = 0;
end;

resid;

model_diagnostics;
steady;
check(qz_zero_threshold=1e-20);

shocks;
var eps_A;
stderr 0.01;
end;

stoch_simul(irf=20,order=1) dY dB SPREAD dCs dCb dCf dHs dHb dq;