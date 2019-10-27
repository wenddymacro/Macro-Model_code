function forecast(var_list,task)

% function forecast(var_list,task)
%   computes mean forecast for a given value of the parameters
%   computes also confidence band for the forecast
%
% INPUTS
%   var_list:    list of variables (character matrix)
%   task:        indicates how to initialize the forecast
%                either 'simul' or 'smoother'
% OUTPUTS
%   nothing is returned but the procedure saves output
%   in oo_.forecast.Mean
%      oo_.forecast.HPDinf
%      oo_.forecast.HPDsup
%
% SPECIAL REQUIREMENTS
%    none
%
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

   global options_ dr_ oo_ M_

   old_options = options_;

   maximum_lag = M_.maximum_lag;
   y_smoothed = oo_.SmoothedVariables;
   options_  = set_default_option(options_,'forecast',12);
   horizon = options_.forecast;
   options_.periods = horizon;
   options_ = set_default_option(options_,'conf_sig',0.9);

   endo_names = M_.endo_names;
   if isempty(var_list)
       var_list = endo_names;
       i_var = 1:M_.endo_nbr;
   else
       i_var = [];
       for i = 1:size(var_list)
           tmp = strmatch(var_list(i,:),endo_names,'exact');
           if isempty(tmp)
               error([var_list(i,:) ' isn''t and endogenous variable'])
           end
           i_var = [i_var; tmp];
       end
   end
   n_var = length(i_var);

   trend = 0;
   switch task
    case 'simul'
     if size(oo_.endo_simul,2) < maximum_lag
         y0 = repmat(oo_.steady_state,1,maximum_lag);
     else
         y0 = oo_.endo_simul(:,1:maximum_lag);
     end
    case 'smoother'
     y0 = zeros(M_.endo_nbr,maximum_lag);
     for i = 1:M_.endo_nbr
         v_name = deblank(M_.endo_names(i,:));
         y0(i,:) = y_smoothed.(v_name)(end-maximum_lag+1:end)+oo_.dr.ys(i);
     end
     gend = size(y0,2);
     if isfield(oo_.Smoother,'TrendCoeffs')
         var_obs = options_.varobs;
         endo_names = M_.endo_names;
         order_var = oo_.dr.order_var;
         i_var_obs = [];
         trend_coeffs = [];
         for i=1:size(var_obs,1)
             tmp = strmatch(var_obs(i,:),endo_names(i_var,:),'exact');
             if ~isempty(tmp)
                 i_var_obs = [ i_var_obs; tmp];
                 trend_coeffs = [trend_coeffs; oo_.Smoother.TrendCoeffs(i)];
             end
         end
         trend = trend_coeffs*(gend+(1-M_.maximum_lag:horizon));
     end
    otherwise
     error('Wrong flag value')
   end

   if M_.exo_det_nbr == 0
       [yf,int_width,int_width2] = forcst(oo_.dr,y0,options_.periods,var_list);
   else
       exo_det_length = size(oo_.exo_det_simul,1);
       if options_.periods > exo_det_length
           ex = zeros(options_.periods,M_.exo_nbr);
           oo_.exo_det_simul = [ oo_.exo_det_simul;...
                   repmat(oo_.exo_det_steady_state',...
                          options_.periods- ...
                          exo_det_length,1)];
           %ex_det_length,1),1)];
       elseif options_.periods < exo_det_length
           ex = zeros(exo_det_length,M_.exo_nbr);
       end
       [yf,int_width] = simultxdet(y0,dr_,ex,oo_.exo_det_simul,...
                                   options_.order,var_list);
   end

   if ~isscalar(trend)
       yf(i_var_obs,:) = yf(i_var_obs,:) + trend;
   end

   for i=1:n_var
       eval(['oo_.forecast.Mean.' var_list(i,:) '= yf(' int2str(i) ',maximum_lag+1:end)'';']);
       eval(['oo_.forecast.HPDinf.' var_list(i,:) '= yf(' int2str(i) ',maximum_lag+1:end)''-' ...
             ' int_width(:,' int2str(i) ');']);
       eval(['oo_.forecast.HPDsup.' var_list(i,:) '= yf(' int2str(i) ',maximum_lag+1:end)''+' ...
             ' int_width(:,' int2str(i) ');']);
       eval(['oo_.forecast.HPDinf2.' var_list(i,:) '= yf(' int2str(i) ',maximum_lag+1:end)''-' ...
             ' int_width2(:,' int2str(i) ');']);
       eval(['oo_.forecast.HPDsup2.' var_list(i,:) '= yf(' int2str(i) ',maximum_lag+1:end)''+' ...
             ' int_width2(:,' int2str(i) ');']);
   end

   for i=1:M_.exo_det_nbr
       eval(['oo_.forecast.Exogenous.' lgx_det_(i,:) '= M_.exo_det_simul(:,' int2str(i) ');']);
   end

   if options_.nograph == 0
       forecast_graphs(var_list);
   end

   options_ = old_options;



