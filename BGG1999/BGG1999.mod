% Dynare Code For
% Bernanke Gertler Gilchrist 
% "The Financial Accelerator in a Quantitative Business-Cycle Framework."
% Nicola Viegi
% Pretoria, 2010

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c i g ce n rk r q k x a h pi rn ;

varexo e_rn e_g e_a;

parameters beta eta alph delt omeg eps G_Y C_Y I_Y Ce_Y Y_N X rho_a sig_a rho_g sig_g psi K_N R gam mu nu thet rho sig kap;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------

C_Y	    =	0.605770191	; // C/Y
Ce_Y	=	0.01	; //Ce/Y
I_Y 	=	0.184229809	; //% I/Y
G_Y 	=	0.2	; //% G/Y
K_N 	=	2.081759973	; //%K/N
Y_N 	=	0.282494996	; //%Y/N
X   	=	1.1	; //%X
beta    =   0.95;
R       =   1/beta;
alph    =   0.2; //capital share of production
eta     =   1/alph;
omeg    =   0.99;
delt    =   0.025;
rho_a   =   0.99;
rho_g   =   0.95;
psi     =   0.25;
Rk      = R + 0.02;
gam     = 1-0.0272;
mu      =   0.12;
thet    =   0.75;
rho     =   0.96;
sig     =   0.20;
kap     =  ((1-thet)/thet)*(1-thet*beta);
eps     =   (1-delt)/((1-delt) + ((alph/X)*(Y_N/K_N)));
nu      =  0.052092347	;

%----------------------------------------------------------------
% 3. Model (the number refers to the equation in the paper)
%----------------------------------------------------------------

model(linear);

%Aggregate Demand

y = C_Y*c + I_Y*i + G_Y*g + Ce_Y*ce;                     //4.14
c = -r + c(+1);                                          //4.15       
ce = n;                                                  //4.16
rk(+1) - r = -nu*(n -(q + k));                           //4.17
rk = (1-eps)*(y - k(-1) - x) + eps*q - q(-1);            //4.18
q = psi*(i - k(-1));                                     //4.19    q(+1) = psi*(i(+1)-k)

%Aggregate Supply 

y = a + alph*k(-1) + (1-alph)*omeg*h;                    //4.20
y - h - x - c = (eta^(-1))*h;                            //4.21
pi = kap*(-x) + beta*pi(+1);                             //4.22

%Evolution of State Variables

k = delt*i + (1-delt)*k(-1);                             //4.23
n = gam*R*K_N*(rk - r(-1)) + r(-1) + n(-1);              //4.24

% Monetary Policy Rules and Shock Preocesses

rn = rho*rn(-1) + sig*pi(-1) - e_rn;                     //4.25 Taylor Rule
g = rho_g*g(-1) + e_g;                                   //4.26
a = rho_a*a(-1) + e_a;                                   //4.27
rn = r + pi(+1);                                         // Fisher Equation
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

check;
steady;

shocks;
var e_g; stderr 0.1;
var e_a; stderr 0.1;
var e_rn; stderr 1.0;
end;

stoch_simul(irf=24);
