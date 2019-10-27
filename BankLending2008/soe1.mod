
var RR_US RR_US_BAR
    UNR_US UNR_US_GAP UNR_US_BAR
    PIE_US PIE_US4 Y_US LGDP_US LGDP_US_BAR RS_US G_US LCPI_US
    E4_PIE_US4 E1_Y_US E1_PIE_US 
    UNR_G_US GROWTH_US GROWTH4_US GROWTH4_US_BAR
    BLT_US BLT_US_BAR
    E E2
;

varexo RES_RR_US_BAR RES_UNR_US_GAP RES_UNR_US_BAR
       RES_RS_US RES_G_US RES_Y_US RES_LGDP_US_BAR RES_PIE_US
       RES_UNR_G_US
       RES_BLT_US RES_BLT_US_BAR
;

parameters rho_us rr_us_bar_ss alpha_us1 alpha_us2
           tau_us growth_us_ss beta_us1 beta_us2 beta_us3 lambda_us1 lambda_us2 gamma_us1 gamma_us2 gamma_us4
           pietar_us_ss 
           alpha_us3
           kappa_us
           theta
;


alpha_us1    = 0.7796;
alpha_us2    = 0.1874;
alpha_us3    = 0.3352;
rho_us       = 0.9;
rr_us_bar_ss = 1.8456;
tau_us      = 0.10;
growth_us_ss= 2.7599;
beta_us1    = 0.9223;
beta_us2    = 0.1198;
beta_us3    = 0.1300;
lambda_us1  = 0.6380;
lambda_us2  = 0.1948;
gamma_us1   = 0.8318;
gamma_us2   = 1.8091;
gamma_us4   = 0.5332;
pietar_us_ss= 2.5;
kappa_us   = 20.0890;
theta = 1;


model;

UNR_US_GAP = alpha_us1*UNR_US_GAP(-1) + alpha_us2*Y_US + RES_UNR_US_GAP;

UNR_US_GAP = UNR_US_BAR - UNR_US;

UNR_US_BAR = UNR_US_BAR(-1) + UNR_G_US + RES_UNR_US_BAR;

UNR_G_US = (1-alpha_us3)*UNR_G_US(-1) + RES_UNR_G_US;




LGDP_US = LGDP_US_BAR + Y_US;

G_US = tau_us*growth_us_ss + (1-tau_us)*G_US(-1) + RES_G_US;

LGDP_US_BAR = LGDP_US_BAR(-1) + G_US/4 + RES_LGDP_US_BAR;

Y_US = beta_us1*Y_US(-1) + beta_us2*Y_US(+1) - beta_us3*(RR_US(-1) - RR_US_BAR(-1)) - 
theta*(0.04*(E(-1)+E(-9))+0.08*(E(-2)+E(-8))+0.12*(E(-3)+E(-7))+0.16*(E(-4)+E(-6))+0.2*E(-5)) + RES_Y_US; 


E = RES_BLT_US;
BLT_US = BLT_US_BAR - kappa_us*Y_US(+4) + RES_BLT_US;
BLT_US_BAR = BLT_US_BAR(-1) + RES_BLT_US_BAR;



GROWTH_US = 4*(LGDP_US - LGDP_US(-1));
GROWTH4_US = LGDP_US - LGDP_US(-4);
GROWTH4_US_BAR = LGDP_US_BAR - LGDP_US_BAR(-4);



PIE_US = lambda_us1*PIE_US4(+4) + (1-lambda_us1)*PIE_US4(-1) + lambda_us2*Y_US(-1) - RES_PIE_US;

RS_US = gamma_us1*RS_US(-1) +
      (1-gamma_us1)*(RR_US_BAR + PIE_US4(+3) + gamma_us2*(PIE_US4(+3)-pietar_us_ss) + gamma_us4*Y_US) + RES_RS_US;


RR_US = RS_US - PIE_US(+1);

LCPI_US = LCPI_US(-1) + PIE_US/4;

RR_US_BAR = rho_us*rr_us_bar_ss + (1-rho_us)*RR_US_BAR(-1) + RES_RR_US_BAR;

PIE_US4 = (PIE_US + PIE_US(-1) + PIE_US(-2) + PIE_US(-3))/4;





// reporting expectations


E4_PIE_US4 = PIE_US4(+4);
E1_PIE_US = PIE_US(+1);
E1_Y_US = Y_US(+1);
E2 = theta*(0.04*(E(-1)+E(-9))+0.08*(E(-2)+E(-8))+0.12*(E(-3)+E(-7))+0.16*(E(-4)+E(-6))+0.2*E(-5));


end;



shocks; 
var RES_RR_US_BAR;  stderr 0.1; 
var RES_UNR_US_GAP; stderr 0.2;
var RES_UNR_US_BAR;  stderr 0.10;
var RES_RS_US	;stderr 0.7;
var RES_G_US	;stderr 0.10;
var RES_Y_US	;stderr 0.25;
var RES_LGDP_US_BAR;stderr 0.05;
var RES_PIE_US    ;stderr 0.7; 
var RES_UNR_G_US ;stderr 0.10;
var RES_BLT_US;    stderr 0.4;
var RES_BLT_US_BAR; stderr 0.2;



//var RES_BLT_US,RES_G_US=(.1*0.40*0.1);
//var RES_BLT_US,RES_LGDP_US_BAR=(.1*0.40*0.05);
var RES_Y_US,RES_G_US=(.1*0.25*0.1);
var RES_LGDP_US_BAR,RES_PIE_US=(.1*0.05*0.7);

end;

unit_root_vars UNR_US_BAR UNR_US LCPI_US LGDP_US LGDP_US_BAR BLT_US BLT_US_BAR;;

steady;

check;


estimated_params;
alpha_us1,    beta_pdf,      0.8, 0.1;
alpha_us2,    gamma_pdf,     0.3, 0.20;
alpha_us3,    beta_pdf,     0.5, 0.20;


growth_us_ss,normal_pdf,    2.5, 0.25;
rr_us_bar_ss, normal_pdf,    2.0, 0.2;

rho_us  ,beta_pdf,      0.9, 0.05;
tau_us  ,beta_pdf,      0.1, 0.05;

beta_us1    ,gamma_pdf,      0.75, 0.10;
beta_us2    ,beta_pdf,      0.15, 0.1; 
beta_us3    ,gamma_pdf,     0.20, 0.05;
lambda_us1  ,beta_pdf,      0.50, 0.10;
lambda_us2  ,gamma_pdf,     0.25, 0.05;
gamma_us1   ,beta_pdf,      0.5, 0.05;
gamma_us2   ,gamma_pdf,     1.50, 0.30;
gamma_us4   ,gamma_pdf,     0.20, 0.05;



//gamma_us1   ,beta_pdf,      0.7, 0.05;
//gamma_us2   ,gamma_pdf,     1.50, 0.30;
//gamma_us4   ,gamma_pdf,     0.50, 0.05;

kappa_us   ,gamma_pdf,     20, 0.5;
theta          ,gamma_pdf,     1, 0.5;



stderr RES_UNR_US_GAP, inv_gamma_pdf, 0.2  , inf;
stderr RES_UNR_US_BAR, inv_gamma_pdf,  0.1  , inf;
stderr RES_UNR_G_US, inv_gamma_pdf,  0.10  ,inf;



stderr RES_Y_US, inv_gamma_pdf,       0.25 , inf; 
stderr RES_LGDP_US_BAR, inv_gamma_pdf, 0.05 , inf; 
stderr RES_G_US, inv_gamma_pdf,  0.10  ,inf;


stderr RES_PIE_US,   inv_gamma_pdf,  0.7  , inf; 
stderr RES_RS_US,	inv_gamma_pdf,  0.7  , inf; 
stderr RES_RR_US_BAR, inv_gamma_pdf,  0.2  , inf;

stderr RES_BLT_US         , inv_gamma_pdf,  0.4  , inf;
stderr RES_BLT_US_BAR     , inv_gamma_pdf,  0.2  , inf;

//corr RES_BLT_US, RES_G_US , beta_pdf, 0.40, 0.15;
//corr RES_BLT_US, RES_LGDP_US_BAR, beta_pdf, 0.4, 0.15;
corr RES_Y_US, RES_G_US, beta_pdf, 0.25, 0.1;
corr RES_LGDP_US_BAR,RES_PIE_US, beta_pdf, 0.05, 0.02;

end;

varobs UNR_US RS_US LCPI_US LGDP_US BLT_US;


observation_trends;
LGDP_US (growth_us_ss/4);
LCPI_US (pietar_us_ss/4);
end;

//% Here is what you need to produce the TeX output:

//% Declaration of two global variables: 
//global M_.exo_names_tex M_.endo_names_tex;

//% Declaration of the tex names of the exogenous variables 
//% (note that the variables must be declared in the same 
//% order than in the global M_.exo_names_):
M_.exo_names_tex = 'RES\_RR\_US\_BAR';
M_.exo_names_tex = strvcat(M_.exo_names_tex, 'RES\_UNR\_US\_GAP');
M_.exo_names_tex = strvcat(M_.exo_names_tex, 'RES\_UNR\_US\_BAR');
M_.exo_names_tex = strvcat(M_.exo_names_tex, 'RES\_RS\_US');
M_.exo_names_tex = strvcat(M_.exo_names_tex, 'RES\_G\_US');
M_.exo_names_tex = strvcat(M_.exo_names_tex, 'RES\_Y\_US');
M_.exo_names_tex = strvcat(M_.exo_names_tex, 'RES\_LGDP\_US\_BAR');
M_.exo_names_tex = strvcat(M_.exo_names_tex, 'RES\_PIE\_US');
M_.exo_names_tex = strvcat(M_.exo_names_tex, 'RES\_UNR\_G\_US');
M_.exo_names_tex = strvcat(M_.exo_names_tex, 'RES\_BLT\_US');
M_.exo_names_tex = strvcat(M_.exo_names_tex, 'RES\_BLT\_US\_BAR');

//% Declaration of the tex names of the endogenous variables 
//% (note that the variables must be declared in the same 
//% order than in the global M_.endo_names_):
M_.endo_names_tex = 'RR\_US';
M_.endo_names_tex = strvcat(M_.endo_names_tex,'RR\_US\_BAR');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'UNR\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'UNR\_US\_GAP');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'UNR\_US\_BAR');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'PIE\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'PIE\_US4');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'Y\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'LGDP\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'LGDP\_US\_BAR');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'RS\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'G\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'LCPI\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'E4\_PIE\_US4');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'E1\_Y\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'E1\_PIE\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'UNR\_G\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'GROWTH\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'GROWTH4\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'GROWTH4\_US\_BAR');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'BLT\_US');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'BLT\_US\_BAR');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'E');
M_.endo_names_tex = strvcat(M_.endo_names_tex,'E2');


//
M_.param_names_tex = '\rho _{us}';
M_.param_names_tex = strvcat(M_.param_names_tex,'rr\_us\_bar\_ss');
M_.param_names_tex = strvcat(M_.param_names_tex,'\alpha _{us1}');
M_.param_names_tex = strvcat(M_.param_names_tex,'\alpha _{us2}');
M_.param_names_tex = strvcat(M_.param_names_tex,'\tau _{us}');
M_.param_names_tex = strvcat(M_.param_names_tex,'');
M_.param_names_tex = strvcat(M_.param_names_tex,'growth\_us\_ss');
M_.param_names_tex = strvcat(M_.param_names_tex,'\beta _{us1}');
M_.param_names_tex = strvcat(M_.param_names_tex,'\beta _{us2}');
M_.param_names_tex = strvcat(M_.param_names_tex,'\beta _{us3}');
M_.param_names_tex = strvcat(M_.param_names_tex,'\lambda _{us1}');
M_.param_names_tex = strvcat(M_.param_names_tex,'\lambda _{us2}');
M_.param_names_tex = strvcat(M_.param_names_tex,'\gamma _{us1}');
M_.param_names_tex = strvcat(M_.param_names_tex,'\gamma _{us2}');
M_.param_names_tex = strvcat(M_.param_names_tex,'\gamma _{us4}');
M_.param_names_tex = strvcat(M_.param_names_tex,'pietar\_us\_ss');
M_.param_names_tex = strvcat(M_.param_names_tex,'\alpha _{us3}');
M_.param_names_tex = strvcat(M_.param_names_tex,'\kappa _{us}');
M_.param_names_tex = strvcat(M_.param_names_tex,'\theta _{us}');


//% Declaration of the tex names of the deep parameters to  
//% be estimated (note that the variables must be declared 
//% in the same order than in the global estim_params_.names):
estim_params_.tex = '\alpha _{us}';
estim_params_.tex = strvcat(estim_params_.tex,'\alpha _{us2}');
estim_params_.tex = strvcat(estim_params_.tex,'\alpha _{us3}');
estim_params_.tex = strvcat(estim_params_.tex,'growth \_us\_ss');
estim_params_.tex = strvcat(estim_params_.tex,'rr\_us\_bar\_ss');
estim_params_.tex = strvcat(estim_params_.tex,'\rho _{us}');
estim_params_.tex = strvcat(estim_params_.tex,'\tau _{us}');
estim_params_.tex = strvcat(estim_params_.tex,'\beta _{us1}');
estim_params_.tex = strvcat(estim_params_.tex,'\beta _{us2}');
estim_params_.tex = strvcat(estim_params_.tex,'\beta _{us3}');
estim_params_.tex = strvcat(estim_params_.tex,'\lambda _{us1}');
estim_params_.tex = strvcat(estim_params_.tex,'\lambda _{us2}');
estim_params_.tex = strvcat(estim_params_.tex,'\gamma _{us1}');
estim_params_.tex = strvcat(estim_params_.tex,'\gamma _{us2}');
estim_params_.tex = strvcat(estim_params_.tex,'\gamma _{us4}');
estim_params_.tex = strvcat(estim_params_.tex,'\kappa _{us}');
estim_params_.tex = strvcat(estim_params_.tex,'\theta _{us}');


estim_params_.tex = strvcat(estim_params_.tex,'RES\_UNR\_US\_GAP');
estim_params_.tex = strvcat(estim_params_.tex,'RES\_UNR\_US\_BAR');
estim_params_.tex = strvcat(estim_params_.tex,'RES\_UNR\_G\_US');
estim_params_.tex = strvcat(estim_params_.tex,'RES\_Y\_US');
estim_params_.tex = strvcat(estim_params_.tex,'RES\_LGDP\_US\_BAR');
estim_params_.tex = strvcat(estim_params_.tex,'RES\_G\_US');
estim_params_.tex = strvcat(estim_params_.tex,'RES\_PIE\_US');
estim_params_.tex = strvcat(estim_params_.tex,'RES\_RS\_US');
estim_params_.tex = strvcat(estim_params_.tex,'RES\_RR\_US\_BAR');
estim_params_.tex = strvcat(estim_params_.tex,'RES\_BLT\_US');
estim_params_.tex = strvcat(estim_params_.tex,'RES\_BLT\_US\_BAR');

estim_params_.tex = strvcat(estim_params_.tex,'CC_RES\_Y\_US\_RES\_G\_US');
estim_params_.tex = strvcat(estim_params_.tex,'CC\_RES\_LGDP\_US\_BAR\_RES\_PIE\_US');



options_.kalman_algo = 5;

//estimation(datafile=data,nobs=56,mode_check, mh_replic=0,mh_jscale=0.35,mh_nblocks=1, filtered_vars,filter_step_ahead=[1:12],tex,forecast=12) Y_US PIE_US4 RS_US UNR_US UNR_US_BAR GROWTH_US GROWTH4_US GROWTH4_US_BAR BLT_US; 

estimation(datafile=data,nobs=56,mode_compute=0,mode_file=soe1_mode, mh_replic=0,mh_jscale=0.35,mh_nblocks=1, filtered_vars,filter_step_ahead=[1:12], forecast=12,conf_sig=0.95) Y_US PIE_US4 RS_US UNR_US UNR_US_BAR GROWTH_US GROWTH4_US GROWTH4_US_BAR BLT_US; 



stoch_simul(order=1) Y_US PIE_US4 RS_US UNR_US UNR_US_BAR GROWTH_US GROWTH4_US GROWTH4_US_BAR BLT_US;


//to make more dynamic for report
data;
nobs = 56;

gl_report
//canada_report_extra__US

%% set firstdate
firstdate = '1994Q1'; 

%% set in-sample ahead forecast steps
fcast_steps = [1 4 8 12];

%% calculate the forecasts
max_step = max(fcast_steps);
[fy,fx,Fy,Fx,Udecomp,ndiffuse] = calc_fcast_all_4a;
//[fy,fx,Fy,Fx,Udecomp,ndiffuse] = calc_fcast_all(max_step);
%% save the forecasts
save([M_.dname '_forecasts.mat'], 'fy', 'fx', 'Fy', 'Fx', 'Udecomp', 'ndiffuse');

disp(sprintf('\n'));disp('Reporting forecast errors');
//report_fcast_errors(fy, Fy, ndiffuse, fcast_steps, firstdate);
report_smoothed_errors_4(fx, Fx, ndiffuse, strvcat('Y_US', 'PIE_US4', 'RS_US', 'UNR_US', 'UNR_US_BAR', 'GROWTH_US', 'GROWTH4_US', 'GROWTH4_US_BAR'), fcast_steps, firstdate);

%% order shocks as they should appear in the decompositions (don't put
%% highly contributing shocks close to each other, since they would get
%% similar colors)
exo_ord = [1 4 2 5 3 6 7 8 9 10];

disp('Reporting dynamic forecasts 1 period ahead');
report_fcast_dynamic_4a(fx, Fx, Udecomp, ndiffuse, 1, 1, exo_ord, firstdate,var_list_);
disp('Reporting dynamic forecasts 4 periods ahead');
report_fcast_dynamic_4a(fx, Fx, Udecomp, ndiffuse, 1, 4, exo_ord, firstdate,var_list_);
disp('Reporting dynamic forecasts 8 periods ahead');
report_fcast_dynamic_4a(fx, Fx, Udecomp, ndiffuse, 1, 8, exo_ord, firstdate,var_list_);
disp('Reporting dynamic forecasts 12 periods ahead');
report_fcast_dynamic_4a(fx, Fx, Udecomp, ndiffuse, 1, 12, exo_ord, firstdate,var_list_);

//forecast([],'smoother');
//forecast_graphs([]);


