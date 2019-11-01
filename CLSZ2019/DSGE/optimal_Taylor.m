% solve for optimal policy rule using fminsearch code
% calculate social welfare in the benchmark case and optimal rules, using second-order approximation 
% it will take  hours to run these codes.
clear

global oo_ M_ options_;

dynare model_Taylor noclearall;
RP1 = -0.397/(1-0.391);
RY1 = 0.183/(1-0.391);
%RY1 = -1.299/(1-0.391);
opt_rule2(RP1,RY1,0,0,0,0)


%% solve for optimal RRR rule
[opt_tau, opt_utility_tau,flag_tau]= fminsearch(@(policy) opt_rule2(RP1,RY1,policy(1),policy(2),0,0),[0 0],optimset('Display','iter', 'TolFun',1e-6))
save respolicy opt_tau opt_utility_tau flag_tau;

[opt_tau1, opt_utility_tau1,flag_tau1]= fminsearch(@(policy) opt_rule2(RP1,RY1,policy(1),policy(2),0,0),[-2 -2],optimset('Display','iter', 'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1;

[opt_tau2, opt_utility_tau2,flag_tau2]= fminsearch(@(policy) opt_rule2(RP1,RY1,policy(1),policy(2),0,0),[2 2],optimset('Display','iter', 'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2;

[opt_tau3, opt_utility_tau3,flag_tau3]= fminsearch(@(policy) opt_rule2(RP1,RY1,policy(1),policy(2),0,0),[5 -5],optimset('Display','iter', 'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3;

[opt_tau4, opt_utility_tau4,flag_tau4]= fminsearch(@(policy) opt_rule2(RP1,RY1,policy(1),policy(2),0,0),[-2 2],optimset('Display','iter', 'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4;

[opt_r, opt_utility_r,flag_r]= fminsearch(@(policy) opt_rule2(policy(1),policy(2),0,0,0,0),[RP1,RY1],optimset('Display','iter', 'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r ;

[opt_r1, opt_utility_r1,flag_r1]= fminsearch(@(policy) opt_rule2(policy(1),policy(2),0,0,0,0),[RP1,-RY1],optimset('Display','iter', 'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1;


[opt_both, opt_utility_both, flag_both]= fminsearch(@(policy) opt_rule2(policy(1),policy(2),policy(3),policy(4),0,0),[RP1 RY1 0 0],optimset('Display','iter', 'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both;

[opt_both1, opt_utility_both1, flag_both1]= fminsearch(@(policy) opt_rule2(policy(1),policy(2),policy(3),policy(4),0,0),[RP1 RY1 opt_tau(1) opt_tau(2)],optimset('Display','iter',  'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both opt_both1 opt_utility_both1 flag_both1;

[opt_both2, opt_utility_both2, flag_both2]= fminsearch(@(policy) opt_rule2(policy(1),policy(2),policy(3),policy(4),0,0),[opt_r(1) opt_r(2) 0 0],optimset('Display','iter', 'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both opt_both1 opt_utility_both1 flag_both1 opt_both2 opt_utility_both2 flag_both2;

[opt_both3, opt_utility_both3, flag_both3]= fminsearch(@(policy) opt_rule2(policy(1),policy(2),policy(3),policy(4),0,0),[opt_r(1) opt_r(2) opt_tau(1) opt_tau(2)],optimset('Display','iter', 'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both opt_both1 opt_utility_both1 flag_both1 opt_both2 opt_utility_both2 flag_both2 opt_both3 opt_utility_both3 flag_both3;

[opt_both4, opt_utility_both4, flag_both4]= fminsearch(@(policy) opt_rule2(policy(1),policy(2),policy(3),policy(4),0,0),[opt_r(1) opt_r(2) opt_tau3(1) opt_tau3(2)],optimset('Display','iter', 'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both opt_both1 opt_utility_both1 flag_both1 opt_both2 opt_utility_both2 flag_both2 opt_both3 opt_utility_both3 flag_both3 opt_both4 opt_utility_both4 flag_both4;

benchmark_utility = opt_rule2(RP1,RY1,0,0,0,0);
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both opt_both1 opt_utility_both1 flag_both1 opt_both2 opt_utility_both2 flag_both2 opt_both3 opt_utility_both3 flag_both3 opt_both4 opt_utility_both4 flag_both4 benchmark_utility;

[opt_both5, opt_utility_both5, flag_both5]= fminsearch(@(policy) opt_rule2(policy(1),policy(2),policy(3),policy(4),0,0),[opt_both4(1) opt_both4(2) opt_both4(3) opt_both4(4)],optimset('Display','iter', 'TolFun',1e-6))
load respolicy
save respolicy opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both opt_both1 opt_utility_both1 flag_both1 opt_both2 opt_utility_both2 flag_both2 opt_both3 opt_utility_both3 flag_both3 opt_both4 opt_utility_both4 flag_both4 benchmark_utility opt_both5 opt_utility_both5 flag_both5;


clear
load shocktype
load respolicy
if token == 1
    save respolicy_Agg opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both opt_both1 opt_utility_both1 flag_both1 opt_both2 opt_utility_both2 flag_both2 opt_both3 opt_utility_both3 flag_both3 opt_both4 opt_utility_both4 flag_both4 benchmark_utility opt_both5 opt_utility_both5 flag_both5;

elseif token == 2
    save respolicy_AS opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both opt_both1 opt_utility_both1 flag_both1 opt_both2 opt_utility_both2 flag_both2 opt_both3 opt_utility_both3 flag_both3 opt_both4 opt_utility_both4 flag_both4 benchmark_utility opt_both5 opt_utility_both5 flag_both5;
[opt_both6, opt_utility_both6, flag_both6]= fminsearch(@(policy) opt_rule2(policy(1),policy(2),policy(3),policy(4),0,0),[-1000 -20 -100 -30],optimset('Display','iter', 'TolFun',1e-6))
load respolicy_AS
save respolicy_AS opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both opt_both1 opt_utility_both1 flag_both1 opt_both2 opt_utility_both2 flag_both2 opt_both3 opt_utility_both3 flag_both3 opt_both4 opt_utility_both4 flag_both4 benchmark_utility opt_both5 opt_utility_both5 flag_both5 opt_both6 opt_utility_both6 flag_both6;
opt_both70 = [  -9.935653908470929  -1.650738014284109   1.532026884655346  -0.309540520982549]*100;% used as initial point to solve opt_both7, very close to opt_both6
[opt_both7, opt_utility_both7, flag_both7]= fminsearch(@(policy) opt_rule2(policy(1),policy(2),policy(3),policy(4),0,0),[opt_both70(1) opt_both70(2) opt_both70(3) opt_both70(4)],optimset('Display','iter', 'TolFun',1e-6,'MaxFunEvals',8000))
load respolicy_AS
save respolicy_AS opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both opt_both1 opt_utility_both1 flag_both1 opt_both2 opt_utility_both2 flag_both2 opt_both3 opt_utility_both3 flag_both3 opt_both4 opt_utility_both4 flag_both4 benchmark_utility opt_both5 opt_utility_both5 flag_both5 opt_both6 opt_utility_both6 flag_both6 opt_both7 opt_utility_both7 flag_both7;

else
    save respolicy_AP opt_tau opt_utility_tau flag_tau opt_tau1 opt_utility_tau1 flag_tau1 opt_tau2 opt_utility_tau2 flag_tau2 opt_tau3 opt_utility_tau3 flag_tau3 opt_tau4 opt_utility_tau4 flag_tau4 opt_r opt_utility_r flag_r opt_r1 opt_utility_r1 flag_r1 opt_both opt_utility_both flag_both opt_both1 opt_utility_both1 flag_both1 opt_both2 opt_utility_both2 flag_both2 opt_both3 opt_utility_both3 flag_both3 opt_both4 opt_utility_both4 flag_both4 benchmark_utility opt_both5 opt_utility_both5 flag_both5;
end
