function out = calcss(x,parameter)

betas =parameter(1);
betab =parameter(2);
jj     =parameter(3);
etaa   =parameter(4);
betaf =parameter(5);
alphaa =parameter(6);
rhoA  =parameter(7);
sigmaA=parameter(8);
kSS     =parameter(9);
CRRSS   =parameter(10);


gammaa = 1-CRRSS;
SIG = (1-betab-(1-betab/betaf*(1-(1-betaf*1/betas)*gammaa))*kSS*((1-(1-betaf*1/betas)*gammaa)*1/betaf)^(-1));
SIG2 = (1/betas*gammaa+1-gammaa-(1-(1-betaf*1/betas)*gammaa)*1/betaf)*kSS*((1-(1-betaf*1/betas)*gammaa)*1/betaf)^(-1);
SIG3 = -gammaa*kSS*jj*(1-1/betas)*(1-alphaa)/alphaa*SIG^(-1)*((1-(1-betaf*1/betas)*gammaa)*1/betaf)^(-1);
NbSS = x;
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

out = [CbSS + RbSS*BSS - BSS - WbSS*NbSS];
