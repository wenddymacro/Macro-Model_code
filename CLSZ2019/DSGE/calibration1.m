function parameter=calibration1
% calibrate MU,XIS,XIP,AP to target ws/wp, ys/yp, SOE bankruptcy ratio, POE
% bankruptcy ratio
clear;
% fixed calibrated values
BETA = 0.995;
ETA = 2;
PSI = 18;
OMEGAK=1; %3;
DELTA = 0.035;
G = 1+0.05/4;
PISS = 1+0.02/4;
GSS = 0.14;


EPS = 10;
OMEGAP=22;

K = 1/0.63;
OMEGAM = (K-1)/K;
ALPHA = 0.5;
AS = 1;
THETA = 0.9;
MS =0.15;%optimal tau increase with M
MP =0.15;

% THETA = 0.95;
% MS = 0.05;
% MP = 0.05;

LS = 1;
LP = 0;
TAU = 0.15;

% target values
%YSP = 0.3;
bcs = 0.25;%optimal tau increase with bc
bcp = 0.1;
% flexible calibrated values 
r = G*PISS/BETA;
x = EPS/(EPS-1);
rs = (r-1)/(1-TAU)+1;
rp = r;
omegas = OMEGAM/((1-bcs/4)^(1/K));
omegap = OMEGAM/((1-bcp/4)^(1/K));
FP = 1/(K-1)*OMEGAM^K*omegap^(-K+1);
GP = LP*omegap+(1-LP)*(1-(1-MP)*K/(K-1))*OMEGAM^K*omegap^(-K+1)+(1-MP)*(1-LP)*K/(K-1)*OMEGAM;
G1P = LP+(1-LP)*(1-MP*K)*OMEGAM^K*omegap^(-K);
F1P = -OMEGAM^K*omegap^(-K);
BP_YP = GP/rp;
NP_YP = -G1P/F1P*FP/rp*PISS*G;


BS_YS = omegas/rs; % given ls = 1
NS_YS = omegas/rs*PISS*G/(K-1); % given ls = 1
FS = 1/(K-1)*OMEGAM^K*omegas^(-K+1);
WE_YP = ALPHA*(1-THETA)*(NP_YP/PISS/G+BP_YP);
WE_YS = ALPHA*(1-THETA)*(NS_YS/PISS/G+BS_YS);


XIP = (NP_YP-WE_YP)/(FP);
XIS = (NS_YS-WE_YS)/(FS);

AP =1.42;

PSIY = 0.45;%optimal tau decrease with PSIY
SIGMAM = 3;

SIGG = 0;
SIGDA = 0;
SIGS=0.0;
SIGAP = 0.0;
SIGAS = 0.0;

load shocktype
if token == 1
    SIGS = 0.01;
elseif token == 2
    SIGAS = 0.01;
else
    SIGAP = 0.01;
end



parameter = [DELTA G OMEGAK PSI PSIY ETA BETA EPS OMEGAP PISS ALPHA THETA XIS XIP K OMEGAM AS AP LS LP TAU GSS SIGMAM MS MP SIGS SIGG SIGAP SIGAS SIGDA];
end