function F = solveff(AS,AP,K,OMEGAM,LS,LP,MS,MP,XIS,XIP,ALPHA,THETA,G,PISS,BETA,TAU,omegas,omegap)
% used by solvess.m to solve for bankruptcy cutoffs omegas omegap using Eq. 16-28.
r = G*PISS/BETA; % Eq. 7
rs = (r-1)/(1-TAU)+1; % Eq. 31
rp = r; % Eq. 34

FS = 1/(K-1)*OMEGAM^K*omegas^(-K+1); % fraction of revenue income going to the SOE borrower f(omegas)
GS = LS*omegas+(1-LS)*(1-(1-MS)*K/(K-1))*OMEGAM^K*omegas^(-K+1)+(1-MS)*(1-LS)*K/(K-1)*OMEGAM; % fraction of revenue income going to the SOE lender g(omegas)
G1S = LS+(1-LS)*(1-MS*K)*OMEGAM^K*omegas^(-K); % first order derivative of g(omegas)
F1S = -OMEGAM^K*omegas^(-K);% first order derivative of f(omegas)
BS_YS = GS/rs; % Eq. 32
NS_YS = -G1S/F1S*FS/rs*PISS*G; % Eq. 27

FP = 1/(K-1)*OMEGAM^K*omegap^(-K+1); % fraction of revenue income going to the POE borrower f(omegap)
GP = LP*omegap+(1-LP)*(1-(1-MP)*K/(K-1))*OMEGAM^K*omegap^(-K+1)+(1-MP)*(1-LP)*K/(K-1)*OMEGAM;% fraction of revenue income going to the SOE lender g(omegap)
G1P = LP+(1-LP)*(1-MP*K)*OMEGAM^K*omegap^(-K);% first order derivative of g(omegap)
F1P = -OMEGAM^K*omegap^(-K);% first order derivative of g(omegap)
BP_YP = GP/rp; % Eq. 26
NP_YP = -G1P/F1P*FP/rp*PISS*G;% Eq. 27


WE_YP = ALPHA*(1-THETA)*(NP_YP/PISS/G+BP_YP);
WE_YS = ALPHA*(1-THETA)*(NS_YS/PISS/G+BS_YS);


F = [WE_YP+XIP*FP-NP_YP;WE_YS+XIS*FS-NS_YS];% return to investment is proportional to TFP (derived from Eq 18, 20-22)  ; Eq 28 for each type of firm 
    

