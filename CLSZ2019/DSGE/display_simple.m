%% use respolicy_X matrix to generate volatilites and welfare results
%addpath d:\dynare\4.4.3\matlab
clear


%% define oo_benchmark oo_r oo_tau oo_both
% oo_benchmark:oo_ from benchmark policy
% oo_r: oo_ file under optimal interest rate policy
% oo_tau: oo_ file under optimal reserve requirement policy
% oo_both: oo_ file under joint optimal policy 
dynare model_Taylor 

load shocktype
if token == 1
load respolicy_Agg % tau4 r0 both4

set_param_value('RP',opt_both4(1));
set_param_value('RY',opt_both4(2));
set_param_value('TP',opt_both4(3));
set_param_value('TY',opt_both4(4));
set_param_value('RR',0);
set_param_value('TT',0);

options_.order =1;
var_list_=[];
info = stoch_simul(var_list_);
oo_both = oo_;
save resdisplay oo_both;

set_param_value('RP',opt_r(1));
set_param_value('RY',opt_r(2));
set_param_value('TP',0);
set_param_value('TY',0);
set_param_value('RR',0);
set_param_value('TT',0);

options_.order =1;
var_list_=[];
info = stoch_simul(var_list_);
oo_r = oo_;
save resdisplay oo_both oo_r;

RP1 = -0.397/(1-0.391);
RY1 = 0.183/(1-0.391);
set_param_value('RP',RP1);
set_param_value('RY',RY1);
set_param_value('TP',opt_tau4(1));
set_param_value('TY',opt_tau4(2));
set_param_value('RR',0);
set_param_value('TT',0);

options_.order =1;
var_list_=[];
info = stoch_simul(var_list_);
oo_tau = oo_;
save resdisplay oo_both oo_r oo_tau;

elseif token == 2
load respolicy_AS % tau2 r1 both7

set_param_value('RP',opt_both7(1));
set_param_value('RY',opt_both7(2));
set_param_value('TP',opt_both7(3));
set_param_value('TY',opt_both7(4));
set_param_value('RR',0);
set_param_value('TT',0);

options_.order =1;
var_list_=[];
info = stoch_simul(var_list_);
oo_both = oo_;
save resdisplay oo_both;

set_param_value('RP',opt_r1(1));
set_param_value('RY',opt_r1(2));
set_param_value('TP',0);
set_param_value('TY',0);
set_param_value('RR',0);
set_param_value('TT',0);

options_.order =1;
var_list_=[];
info = stoch_simul(var_list_);
oo_r = oo_;
save resdisplay oo_both oo_r;

RP1 = -0.397/(1-0.391);
RY1 = 0.183/(1-0.391);
set_param_value('RP',RP1);
set_param_value('RY',RY1);
set_param_value('TP',opt_tau2(1));
set_param_value('TY',opt_tau2(2));
set_param_value('RR',0);
set_param_value('TT',0);

options_.order =1;
var_list_=[];
info = stoch_simul(var_list_);
oo_tau = oo_;
save resdisplay oo_both oo_r oo_tau;

    
else
load respolicy_AP % tau3 r1 both3

set_param_value('RP',opt_both3(1));
set_param_value('RY',opt_both3(2));
set_param_value('TP',opt_both3(3));
set_param_value('TY',opt_both3(4));
set_param_value('RR',0);
set_param_value('TT',0);

options_.order =1;
var_list_=[];
info = stoch_simul(var_list_);
oo_both = oo_;
save resdisplay oo_both;

set_param_value('RP',opt_r1(1));
set_param_value('RY',opt_r1(2));
set_param_value('TP',0);
set_param_value('TY',0);
set_param_value('RR',0);
set_param_value('TT',0);

options_.order =1;
var_list_=[];
info = stoch_simul(var_list_);
oo_r = oo_;
save resdisplay oo_both oo_r;

RP1 = -0.397/(1-0.391);
RY1 = 0.183/(1-0.391);
set_param_value('RP',RP1);
set_param_value('RY',RY1);
set_param_value('TP',opt_tau3(1));
set_param_value('TY',opt_tau3(2));
set_param_value('RR',0);
set_param_value('TT',0);

options_.order =1;
var_list_=[];
info = stoch_simul(var_list_);
oo_tau = oo_;
save resdisplay oo_both oo_r oo_tau;

end



RP1 = -0.397/(1-0.391);
RY1 = 0.183/(1-0.391);
set_param_value('RP',RP1);
set_param_value('RY',RY1);
set_param_value('TP',0);
set_param_value('TY',0);
set_param_value('RR',0);
set_param_value('TT',0);
options_.order =1;
var_list_=[];
info = stoch_simul(var_list_);
oo_benchmark = oo_;
save resdisplay oo_both oo_r oo_tau oo_benchmark;

%% obtain steady state mean and variance from oo_ file
clear
load shocktype
if token == 1
    load respolicy_Agg;
elseif token == 2
    load respolicy_AS;
else
    load respolicy_AP;
end

load resdisplay
var_benchmark_k=oo_benchmark.var(1,1);
var_benchmark_i=oo_benchmark.var(2,2);
var_benchmark_uc=oo_benchmark.var(3,3);
var_benchmark_c=oo_benchmark.var(4,4);
var_benchmark_h=oo_benchmark.var(5,5);
var_benchmark_w=oo_benchmark.var(6,6);
var_benchmark_r=oo_benchmark.var(7,7);
var_benchmark_pi=oo_benchmark.var(8,8);
var_benchmark_uk=oo_benchmark.var(9,9);
var_benchmark_rk=oo_benchmark.var(10,10);
var_benchmark_x=oo_benchmark.var(11,11);
var_benchmark_y=oo_benchmark.var(12,12);
var_benchmark_wes=oo_benchmark.var(13,13);
var_benchmark_wep=oo_benchmark.var(14,14);
var_benchmark_ys=oo_benchmark.var(15,15);
var_benchmark_ks=oo_benchmark.var(16,16);
var_benchmark_hes=oo_benchmark.var(17,17);
var_benchmark_hhs=oo_benchmark.var(18,18);
var_benchmark_ns=oo_benchmark.var(19,19);
var_benchmark_bs=oo_benchmark.var(20,20);
var_benchmark_as=oo_benchmark.var(21,21);
var_benchmark_rs=oo_benchmark.var(22,22);
var_benchmark_omegas=oo_benchmark.var(23,23);
var_benchmark_yp=oo_benchmark.var(24,24);
var_benchmark_kp=oo_benchmark.var(25,25);
var_benchmark_hep=oo_benchmark.var(26,26);
var_benchmark_hhp=oo_benchmark.var(27,27);
var_benchmark_np=oo_benchmark.var(28,28);
var_benchmark_bp=oo_benchmark.var(29,29);
var_benchmark_ap=oo_benchmark.var(30,30);
var_benchmark_rp=oo_benchmark.var(31,31);
var_benchmark_omegap=oo_benchmark.var(32,32);
var_benchmark_z_s=oo_benchmark.var(33,33);
var_benchmark_tau=oo_benchmark.var(34,34);
var_benchmark_gdp=oo_benchmark.var(35,35);
var_benchmark_ytfp=oo_benchmark.var(36,36);
var_benchmark_gtfp=oo_benchmark.var(37,37);
var_benchmark_rzs=oo_benchmark.var(38,38);
var_benchmark_rzp=oo_benchmark.var(39,39);
var_benchmark_utility=oo_benchmark.var(40,40);
var_benchmark_totutility=oo_benchmark.var(41,41);
var_benchmark_g=oo_benchmark.var(42,42);
var_benchmark_bcs=oo_benchmark.var(43,43);
var_benchmark_bcp=oo_benchmark.var(44,44);
var_benchmark_q=oo_benchmark.var(45,45);
var_benchmark_levs=oo_benchmark.var(46,46);
var_benchmark_levp=oo_benchmark.var(47,47);
var_benchmark_omegas1=oo_benchmark.var(48,48);
var_benchmark_omegap1=oo_benchmark.var(49,49);
var_benchmark_ps=oo_benchmark.var(50,50);
var_benchmark_pp=oo_benchmark.var(51,51);

var_tau_k=oo_tau.var(1,1);
var_tau_i=oo_tau.var(2,2);
var_tau_uc=oo_tau.var(3,3);
var_tau_c=oo_tau.var(4,4);
var_tau_h=oo_tau.var(5,5);
var_tau_w=oo_tau.var(6,6);
var_tau_r=oo_tau.var(7,7);
var_tau_pi=oo_tau.var(8,8);
var_tau_uk=oo_tau.var(9,9);
var_tau_rk=oo_tau.var(10,10);
var_tau_x=oo_tau.var(11,11);
var_tau_y=oo_tau.var(12,12);
var_tau_wes=oo_tau.var(13,13);
var_tau_wep=oo_tau.var(14,14);
var_tau_ys=oo_tau.var(15,15);
var_tau_ks=oo_tau.var(16,16);
var_tau_hes=oo_tau.var(17,17);
var_tau_hhs=oo_tau.var(18,18);
var_tau_ns=oo_tau.var(19,19);
var_tau_bs=oo_tau.var(20,20);
var_tau_as=oo_tau.var(21,21);
var_tau_rs=oo_tau.var(22,22);
var_tau_omegas=oo_tau.var(23,23);
var_tau_yp=oo_tau.var(24,24);
var_tau_kp=oo_tau.var(25,25);
var_tau_hep=oo_tau.var(26,26);
var_tau_hhp=oo_tau.var(27,27);
var_tau_np=oo_tau.var(28,28);
var_tau_bp=oo_tau.var(29,29);
var_tau_ap=oo_tau.var(30,30);
var_tau_rp=oo_tau.var(31,31);
var_tau_omegap=oo_tau.var(32,32);
var_tau_z_s=oo_tau.var(33,33);
var_tau_tau=oo_tau.var(34,34);
var_tau_gdp=oo_tau.var(35,35);
var_tau_ytfp=oo_tau.var(36,36);
var_tau_gtfp=oo_tau.var(37,37);
var_tau_rzs=oo_tau.var(38,38);
var_tau_rzp=oo_tau.var(39,39);
var_tau_utility=oo_tau.var(40,40);
var_tau_totutility=oo_tau.var(41,41);
var_tau_g=oo_tau.var(42,42);
var_tau_bcs=oo_tau.var(43,43);
var_tau_bcp=oo_tau.var(44,44);
var_tau_q=oo_tau.var(45,45);
var_tau_levs=oo_tau.var(46,46);
var_tau_levp=oo_tau.var(47,47);
var_tau_omegas1=oo_tau.var(48,48);
var_tau_omegap1=oo_tau.var(49,49);
var_tau_ps=oo_tau.var(50,50);
var_tau_pp=oo_tau.var(51,51);

var_r_k=oo_r.var(1,1);
var_r_i=oo_r.var(2,2);
var_r_uc=oo_r.var(3,3);
var_r_c=oo_r.var(4,4);
var_r_h=oo_r.var(5,5);
var_r_w=oo_r.var(6,6);
var_r_r=oo_r.var(7,7);
var_r_pi=oo_r.var(8,8);
var_r_uk=oo_r.var(9,9);
var_r_rk=oo_r.var(10,10);
var_r_x=oo_r.var(11,11);
var_r_y=oo_r.var(12,12);
var_r_wes=oo_r.var(13,13);
var_r_wep=oo_r.var(14,14);
var_r_ys=oo_r.var(15,15);
var_r_ks=oo_r.var(16,16);
var_r_hes=oo_r.var(17,17);
var_r_hhs=oo_r.var(18,18);
var_r_ns=oo_r.var(19,19);
var_r_bs=oo_r.var(20,20);
var_r_as=oo_r.var(21,21);
var_r_rs=oo_r.var(22,22);
var_r_omegas=oo_r.var(23,23);
var_r_yp=oo_r.var(24,24);
var_r_kp=oo_r.var(25,25);
var_r_hep=oo_r.var(26,26);
var_r_hhp=oo_r.var(27,27);
var_r_np=oo_r.var(28,28);
var_r_bp=oo_r.var(29,29);
var_r_ap=oo_r.var(30,30);
var_r_rp=oo_r.var(31,31);
var_r_omegap=oo_r.var(32,32);
var_r_z_s=oo_r.var(33,33);
var_r_tau=oo_r.var(34,34);
var_r_gdp=oo_r.var(35,35);
var_r_ytfp=oo_r.var(36,36);
var_r_gtfp=oo_r.var(37,37);
var_r_rzs=oo_r.var(38,38);
var_r_rzp=oo_r.var(39,39);
var_r_utility=oo_r.var(40,40);
var_r_totutility=oo_r.var(41,41);
var_r_g=oo_r.var(42,42);
var_r_bcs=oo_r.var(43,43);
var_r_bcp=oo_r.var(44,44);
var_r_q=oo_r.var(45,45);
var_r_levs=oo_r.var(46,46);
var_r_levp=oo_r.var(47,47);
var_r_omegas1=oo_r.var(48,48);
var_r_omegap1=oo_r.var(49,49);
var_r_ps=oo_r.var(50,50);
var_r_pp=oo_r.var(51,51);

var_both_k=oo_both.var(1,1);
var_both_i=oo_both.var(2,2);
var_both_uc=oo_both.var(3,3);
var_both_c=oo_both.var(4,4);
var_both_h=oo_both.var(5,5);
var_both_w=oo_both.var(6,6);
var_both_r=oo_both.var(7,7);
var_both_pi=oo_both.var(8,8);
var_both_uk=oo_both.var(9,9);
var_both_rk=oo_both.var(10,10);
var_both_x=oo_both.var(11,11);
var_both_y=oo_both.var(12,12);
var_both_wes=oo_both.var(13,13);
var_both_wep=oo_both.var(14,14);
var_both_ys=oo_both.var(15,15);
var_both_ks=oo_both.var(16,16);
var_both_hes=oo_both.var(17,17);
var_both_hhs=oo_both.var(18,18);
var_both_ns=oo_both.var(19,19);
var_both_bs=oo_both.var(20,20);
var_both_as=oo_both.var(21,21);
var_both_rs=oo_both.var(22,22);
var_both_omegas=oo_both.var(23,23);
var_both_yp=oo_both.var(24,24);
var_both_kp=oo_both.var(25,25);
var_both_hep=oo_both.var(26,26);
var_both_hhp=oo_both.var(27,27);
var_both_np=oo_both.var(28,28);
var_both_bp=oo_both.var(29,29);
var_both_ap=oo_both.var(30,30);
var_both_rp=oo_both.var(31,31);
var_both_omegap=oo_both.var(32,32);
var_both_z_s=oo_both.var(33,33);
var_both_tau=oo_both.var(34,34);
var_both_gdp=oo_both.var(35,35);
var_both_ytfp=oo_both.var(36,36);
var_both_gtfp=oo_both.var(37,37);
var_both_rzs=oo_both.var(38,38);
var_both_rzp=oo_both.var(39,39);
var_both_utility=oo_both.var(40,40);
var_both_totutility=oo_both.var(41,41);
var_both_g=oo_both.var(42,42);
var_both_bcs=oo_both.var(43,43);
var_both_bcp=oo_both.var(44,44);
var_both_q=oo_both.var(45,45);
var_both_levs=oo_both.var(46,46);
var_both_levp=oo_both.var(47,47);
var_both_omegas1=oo_both.var(48,48);
var_both_omegap1=oo_both.var(49,49);
var_both_ps=oo_both.var(50,50);
var_both_pp=oo_both.var(51,51);

%% generate table 
BETA = 0.995;
std_gdp =[sqrt(var_benchmark_gdp) sqrt(var_tau_gdp) sqrt(var_r_gdp) sqrt(var_both_gdp)];
std_pi =[sqrt(var_benchmark_pi) sqrt(var_tau_pi) sqrt(var_r_pi) sqrt(var_both_pi)];
std_c =[sqrt(var_benchmark_c) sqrt(var_tau_c) sqrt(var_r_c) sqrt(var_both_c)];
std_h =[sqrt(var_benchmark_h) sqrt(var_tau_h) sqrt(var_r_h) sqrt(var_both_h)];
std_r =[sqrt(var_benchmark_r) sqrt(var_tau_r) sqrt(var_r_r) sqrt(var_both_r)];
std_ys =[sqrt(var_benchmark_ys) sqrt(var_tau_ys) sqrt(var_r_ys) sqrt(var_both_ys)];
std_yp =[sqrt(var_benchmark_yp) sqrt(var_tau_yp) sqrt(var_r_yp) sqrt(var_both_yp)];

RP1 = -0.397/(1-0.391);
RY1 = 0.183/(1-0.391);
load shocktype
if token == 1 % tau4 r0 both4
utility = (exp([0 -opt_utility_tau4-(-benchmark_utility) -opt_utility_r-(-benchmark_utility) -opt_utility_both4-(-benchmark_utility)]*(1-BETA))-1);
format short
policy_table = [RP1 RP1 opt_r(1) opt_both4(1); RY1 RY1 opt_r(2) opt_both4(2);0 opt_tau4(1) 0 opt_both4(3);0 opt_tau4(2) 0 opt_both4(4)]
format long
vol_table = [std_gdp;std_pi;std_c;std_ys;std_yp;utility]

elseif token == 2  % tau2 r1 both6
utility = (exp([0 -opt_utility_tau2-(-benchmark_utility) -opt_utility_r1-(-benchmark_utility) -opt_utility_both7-(-benchmark_utility)]*(1-BETA))-1);
format short
policy_table = [RP1 RP1 opt_r1(1) opt_both7(1); RY1 RY1 opt_r1(2) opt_both7(2);0 opt_tau2(1) 0 opt_both7(3);0 opt_tau2(2) 0 opt_both7(4)]
format long
vol_table = [std_gdp;std_pi;std_c;std_h;std_r;std_ys;std_yp;utility]

else %% tau3 r1 both3
utility = (exp([0 -opt_utility_tau3-(-benchmark_utility) -opt_utility_r1-(-benchmark_utility) -opt_utility_both3-(-benchmark_utility)]*(1-BETA))-1);
format short
policy_table = [RP1 RP1 opt_r1(1) opt_both3(1); RY1 RY1 opt_r1(2) opt_both3(2);0 opt_tau3(1) 0 opt_both3(3);0 opt_tau3(2) 0 opt_both3(4)]
format long
vol_table = [std_gdp;std_pi;std_c;std_ys;std_yp;utility]

end

