function welfare = opt_rule2(RPv,RYv,TPv,TYv,RRv,TTv)
% this m file sovles for social welfare given policy parameters 
% BPHv: mb response to housing prices
% HPHv: mh repsonse to housing prices

global M_ oo_ options_ 

try 
opt_rule1;

var_list_=[];
var_list_ = 'totutility';
var_list_ = char(var_list_, 'tau');

set_param_value('RP',RPv);
set_param_value('RY',RYv);
set_param_value('TP',TPv);
set_param_value('TY',TYv);
set_param_value('RR',RRv);
set_param_value('TT',TTv);

set_dynare_seed('default');
options_.order = 2;
info = stoch_simul(var_list_);
%uncondition welfare
welfare =  -oo_.mean(1);
%condition welfare
%W_pos=strmatch('Wel',M_.endo_names,'exact');
%wlfare=oo_.dr.ys(W_pos)+0.5*oo_.dr.ghs2(oo_.dr.inv_order_var(W_pos))

load result1
if length(find(abs(real(oo_.dr.eigval))>1+1e-3))~=7   | oo_.mean(2)>1 | oo_.mean(2)<0 | resirf == 0 | abs(RPv)>1000 |abs(RYv)>1000 | abs(TPv)>1000| abs(TYv)>1000
 welfare = 1e10;
end;
% | length(find(abs(real(oo_.dr.eigval))<1-2e-3))~=9|max(abs(imag(oo_.dr.eigval)))>1e-6
% length(find(abs(real(oo_.dr.eigval))>1+1e-3))~=7 
catch 
    welfare = 1e10;
end
%| max(abs([RPv RYv TPv TYv]))>100