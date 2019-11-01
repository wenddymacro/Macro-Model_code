global M_ oo_ options_ 

set_param_value('RP',RPv);
set_param_value('RY',RYv);
set_param_value('TP',TPv);
set_param_value('TY',TYv);
set_param_value('RR',RRv);
set_param_value('TT',TTv);

options_.order = 1;
var_list_=[];
info = stoch_simul(var_list_);



parameter=calibration1;
SIGS=parameter(26);
SIGG=parameter(27);
SIGAP=parameter(28);
SIGAS=parameter(29);
SIGDA=parameter(30);

TAU=parameter(21);


resirf = 1;
if SIGAP>0 &(max(-oo_.irfs.tau_eps_ap+TAU)>1|max(oo_.irfs.tau_eps_ap+TAU)>1|min(-oo_.irfs.tau_eps_ap+TAU)<0|min(oo_.irfs.tau_eps_ap+TAU)<0)
resirf = 0;
end;
if SIGAS>0 &(max(-oo_.irfs.tau_eps_as+TAU)>1|max(oo_.irfs.tau_eps_as+TAU)>1|min(-oo_.irfs.tau_eps_as+TAU)<0|min(oo_.irfs.tau_eps_as+TAU)<0)
resirf = 0;
end;
if SIGDA>0 &(max(-oo_.irfs.tau_eps_da+TAU)>1|max(oo_.irfs.tau_eps_da+TAU)>1|min(-oo_.irfs.tau_eps_da+TAU)<0|min(oo_.irfs.tau_eps_da+TAU)<0)
resirf = 0;
end;
if SIGS>0 &(max(-oo_.irfs.tau_eps_s+TAU)>1|max(oo_.irfs.tau_eps_s+TAU)>1|min(-oo_.irfs.tau_eps_s+TAU)<0|min(oo_.irfs.tau_eps_s+TAU)<0)
resirf = 0;
end;
if SIGG>0 &(max(-oo_.irfs.tau_eps_g+TAU)>1|max(oo_.irfs.tau_eps_g+TAU)>1|min(-oo_.irfs.tau_eps_g+TAU)<0|min(oo_.irfs.tau_eps_g+TAU)<0)
resirf = 0;
end;

save result1 resirf;

