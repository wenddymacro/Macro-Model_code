%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

if isoctave || matlab_ver_less_than('8.6')
    clear all
else
    clearvars -global
    clear_persistent_variables(fileparts(which('dynare')), false)
end
tic0 = tic;
% Save empty dates and dseries objects in memory.
dates('initialize');
dseries('initialize');
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info ys0_ ex0_
options_ = [];
M_.fname = 'augmentedrep';
M_.dynare_version = '4.5.5';
oo_.dynare_version = '4.5.5';
options_.dynare_version = '4.5.5';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('augmentedrep.log');
M_.exo_names = 'e_R';
M_.exo_names_tex = 'e\_R';
M_.exo_names_long = 'e_R';
M_.exo_names = char(M_.exo_names, 'e_g');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_g');
M_.exo_names_long = char(M_.exo_names_long, 'e_g');
M_.exo_names = char(M_.exo_names, 'e_z');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_z');
M_.exo_names_long = char(M_.exo_names_long, 'e_z');
M_.exo_names = char(M_.exo_names, 'e_v');
M_.exo_names_tex = char(M_.exo_names_tex, 'e\_v');
M_.exo_names_long = char(M_.exo_names_long, 'e_v');
M_.endo_names = 'x';
M_.endo_names_tex = 'x';
M_.endo_names_long = 'x';
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names_long = char(M_.endo_names_long, 'R');
M_.endo_names = char(M_.endo_names, 'pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi');
M_.endo_names_long = char(M_.endo_names_long, 'pi');
M_.endo_names = char(M_.endo_names, 'g');
M_.endo_names_tex = char(M_.endo_names_tex, 'g');
M_.endo_names_long = char(M_.endo_names_long, 'g');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names_long = char(M_.endo_names_long, 'z');
M_.endo_names = char(M_.endo_names, 'xobs');
M_.endo_names_tex = char(M_.endo_names_tex, 'xobs');
M_.endo_names_long = char(M_.endo_names_long, 'xobs');
M_.endo_names = char(M_.endo_names, 'Robs');
M_.endo_names_tex = char(M_.endo_names_tex, 'Robs');
M_.endo_names_long = char(M_.endo_names_long, 'Robs');
M_.endo_names = char(M_.endo_names, 'pieobs');
M_.endo_names_tex = char(M_.endo_names_tex, 'pieobs');
M_.endo_names_long = char(M_.endo_names_long, 'pieobs');
M_.endo_names = char(M_.endo_names, 'Epi');
M_.endo_names_tex = char(M_.endo_names_tex, 'Epi');
M_.endo_names_long = char(M_.endo_names_long, 'Epi');
M_.endo_names = char(M_.endo_names, 'w');
M_.endo_names_tex = char(M_.endo_names_tex, 'w');
M_.endo_names_long = char(M_.endo_names_long, 'w');
M_.endo_partitions = struct();
M_.param_names = 'Tau';
M_.param_names_tex = 'Tau';
M_.param_names_long = 'Tau';
M_.param_names = char(M_.param_names, 'kappa');
M_.param_names_tex = char(M_.param_names_tex, 'kappa');
M_.param_names_long = char(M_.param_names_long, 'kappa');
M_.param_names = char(M_.param_names, 'rho_R');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_R');
M_.param_names_long = char(M_.param_names_long, 'rho_R');
M_.param_names = char(M_.param_names, 'rho_g');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_g');
M_.param_names_long = char(M_.param_names_long, 'rho_g');
M_.param_names = char(M_.param_names, 'rho_z');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_z');
M_.param_names_long = char(M_.param_names_long, 'rho_z');
M_.param_names = char(M_.param_names, 'psi1');
M_.param_names_tex = char(M_.param_names_tex, 'psi1');
M_.param_names_long = char(M_.param_names_long, 'psi1');
M_.param_names = char(M_.param_names, 'psi2');
M_.param_names_tex = char(M_.param_names_tex, 'psi2');
M_.param_names_long = char(M_.param_names_long, 'psi2');
M_.param_names = char(M_.param_names, 'rstar');
M_.param_names_tex = char(M_.param_names_tex, 'rstar');
M_.param_names_long = char(M_.param_names_long, 'rstar');
M_.param_names = char(M_.param_names, 'piestar');
M_.param_names_tex = char(M_.param_names_tex, 'piestar');
M_.param_names_long = char(M_.param_names_long, 'piestar');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 4;
M_.endo_nbr = 10;
M_.param_nbr = 10;
M_.orig_endo_nbr = 10;
M_.aux_vars = [];
options_.varobs = cell(1);
options_.varobs(1)  = {'xobs'};
options_.varobs(2)  = {'Robs'};
options_.varobs(3)  = {'pieobs'};
options_.varobs_id = [ 6 7 8  ];
M_.Sigma_e = zeros(4, 4);
M_.Correlation_matrix = eye(4, 4);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
M_.det_shocks = [];
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
M_.hessian_eq_zero = 1;
erase_compiled_function('augmentedrep_static');
erase_compiled_function('augmentedrep_dynamic');
M_.orig_eq_nbr = 10;
M_.eq_nbr = 10;
M_.ramsey_eq_nbr = 0;
M_.set_auxiliary_variables = exist(['./' M_.fname '_set_auxiliary_variables.m'], 'file') == 2;
M_.lead_lag_incidence = [
 0 6 16;
 1 7 0;
 0 8 17;
 2 9 0;
 3 10 0;
 0 11 0;
 0 12 0;
 0 13 0;
 4 14 0;
 5 15 0;]';
M_.nstatic = 3;
M_.nfwrd   = 2;
M_.npred   = 5;
M_.nboth   = 0;
M_.nsfwrd   = 2;
M_.nspred   = 5;
M_.ndynamic   = 7;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:4];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(4, 1);
M_.params = NaN(10, 1);
M_.NNZDerivatives = [34; 0; -1];
close all;
estim_params_.var_exo = [];
estim_params_.var_endo = [];
estim_params_.corrx = [];
estim_params_.corrn = [];
estim_params_.param_vals = [];
estim_params_.param_vals = [estim_params_.param_vals; 10, 1.2, (-Inf), Inf, 5, NaN, NaN, 0, 2, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 6, 0.8, NaN, NaN, 2, 1.1, 0.5, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 7, NaN, NaN, NaN, 2, 0.25, 0.15, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 3, NaN, NaN, NaN, 1, 0.5, 0.2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 9, NaN, NaN, NaN, 2, 4, 2, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 8, NaN, NaN, NaN, 2, 2, 1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 2, NaN, NaN, NaN, 2, 0.5, 0.35, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 1, NaN, NaN, NaN, 2, 2, 0.5, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 4, NaN, NaN, NaN, 1, 0.7, 0.1, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 5, NaN, NaN, NaN, 1, 0.7, 0.1, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 1, NaN, NaN, NaN, 4, 0.31, 0.16, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 2, NaN, NaN, NaN, 4, 0.38, 0.2, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 3, NaN, NaN, NaN, 4, 1, 0.52, NaN, NaN, NaN ];
estim_params_.corrx = [estim_params_.corrx; 3 2, NaN, (-1), 1, 3, 0, 0.4, NaN, NaN, NaN ];
estim_params_.corrx = [estim_params_.corrx; 1 2, NaN, (-1), 1, 3, 0, 0.4, NaN, NaN, NaN ];
estim_params_.corrx = [estim_params_.corrx; 3 1, NaN, (-1), 1, 3, 0, 0.4, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 4, NaN, NaN, NaN, 4, 2, 0.5, NaN, NaN, NaN ];
estim_params_.corrx = [estim_params_.corrx; 1 4, NaN, (-1), 1, 3, 0, 0.1, NaN, NaN, NaN ];
estim_params_.corrx = [estim_params_.corrx; 2 4, NaN, (-1), 1, 3, 0, 0.1, NaN, NaN, NaN ];
estim_params_.corrx = [estim_params_.corrx; 3 4, NaN, (-1), 1, 3, (-0.7), 0.1, NaN, NaN, NaN ];
options_.mh_drop = 0.5;
options_.mh_jscale = 0.4;
options_.mh_replic = 250000;
options_.mode_check.status = 1;
options_.mode_compute = 6;
options_.datafile = 'nicoloindeter';
options_.xls_range = 'A1:C501';
options_.graph_format = char('pdf','fig');
options_.order = 1;
var_list_ = char();
oo_recursive_=dynare_estimation(var_list_);
save('augmentedrep_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('augmentedrep_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('augmentedrep_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('augmentedrep_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('augmentedrep_results.mat', 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save('augmentedrep_results.mat', 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save('augmentedrep_results.mat', 'oo_recursive_', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc(tic0)) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
