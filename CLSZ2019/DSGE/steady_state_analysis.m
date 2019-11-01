clear
% effect of required reserve ratio on the steady state welfare
token = 1;
save shocktype token;
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


vTAU = [0.01:0.01:0.5];
for j = 1:1:length(vTAU)
options = optimset('MaxFunEvals',100000,'MaxIter',10000,'TolFun',1e-16);
[FF,fval,exitflag(j)] =fsolve(@(x)solveff(AS,AP,K,OMEGAM,LS,LP,MS,MP,XIS,XIP,ALPHA,THETA,G,PISS,BETA,vTAU(j),OMEGAM*(1+exp(x(1))),OMEGAM*(1+exp(x(2)))),[log(1.02),log(1.01)],options);
[k(j), i(j), uc(j)      , c(j), h(j) ,w(j),r(j) ,pi(j), uk(j), rk(j), x(j), y(j), wes(j),wep(j), ys(j), ks(j), hes(j), hhs(j), ns(j), bs(j), as(j) ,rs(j), omegas(j), yp(j), kp(j), hep(j), hhp(j), np(j), bp(j), ap(j), rp(j), omegap(j), tau(j),gdp(j),ytfp(j),gtfp(j),rzs(j),rzp(j),g(j),bcs(j),bcp(j),ps(j),pp(j)] = solvesslog(DELTA, G, OMEGAK, PSI, ETA, BETA, EPS, OMEGAP, PISS, ALPHA, THETA, MS,MP, XIS,XIP, K, OMEGAM, AS,AP, LS, LP, vTAU(j),GSS,SIGMAM,PSIY);
lp(j)=bp(j)-np(j);%POEs' leverage ratio
ls(j)=bs(j)-ns(j);%SOEs' leverage ratio
utility(j) = (1/(1-BETA))*(c(j)-PSI*exp(h(j))^(1+ETA)/(1+ETA));
end
[M,I] = max(utility)
opt_tau = (I-1)/100;

ysp = ys+ps-yp-pp;% SOE to POE output ratio


figure;
subplot(2,2,1)
plot(vTAU(1:50),exp(ysp(1:50)),'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
title('SOE output/POE output','FontName','Times New Roman');
subplot(2,2,2)
plot(vTAU(1:50),exp(ytfp(1:50)),'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
title('Output-based TFP','FontName','Times New Roman');
subplot(2,2,3)
plot(vTAU(1:50),exp(bcs(1:50)),'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
title('SOE bankruptcy ratio','FontName','Times New Roman');
subplot(2,2,4)
plot(vTAU(1:50),100*(exp((1-BETA)*(utility(1:50)-utility(15)))-1),'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
title('Welfare gains (%)','FontName','Times New Roman');

print -dpdf outfig_ss.pdf

