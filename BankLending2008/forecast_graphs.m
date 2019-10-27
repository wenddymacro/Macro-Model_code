function forecast_graphs(var_list)

   global M_ oo_ options_

   nc = 3;
   nr = 3;
   x = 2007 + [0:options_.forecast+3]'/4;
   exo_nbr = M_.exo_nbr;
   endo_names = M_.endo_names;
   fname = M_.fname;
% $$$     varobs = options_.varobs;
    y = oo_.SmoothedVariables;
    ys = oo_.dr.ys;
% $$$     gend = size(y,2);
   yf = oo_.forecast.Mean;
   hpdinf = oo_.forecast.HPDinf;
   hpdsup = oo_.forecast.HPDsup;
   hpdinf2 = oo_.forecast.HPDinf2;
   hpdsup2 = oo_.forecast.HPDsup2;
   if isempty(var_list)
       varlist = endo_names;
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
   nvar = length(i_var);

% $$$     % build trend for smoothed variables if necessary
% $$$     trend = zeros(size(varobs,1),10);
% $$$     if isfield(oo_.Smoother,'TrendCoeffs')
% $$$         trend_coeffs = oo_.Smoother.TrendCoeffs;
% $$$         trend = trend_coeffs*(gend-9:gend);
% $$$     end

   % create subdirectory <fname>/graphs if id doesn't exist
   if ~exist([fname '/graphs'])
       mkdir(fname,'graphs');
   end

   m = 1;
   n_fig = 1;
   figure('Name','Forecasts (I)')
   for j= 1:nvar
       if m > nc*nr;
           eval(['print -depsc ' fname '/graphs/forcst' int2str(n_fig) ...
                '.eps;'])
           n_fig =n_fig+1;
           eval(['figure(''Name'',''Forecast (' int2str(n_fig) ')'');']);
           m = 1;
       end
       subplot(nr,nc,m);
       vn = deblank(endo_names(i_var(j),:));
       obs = 4;
% $$$         k = strmatch(vn,varobs,'exact');
% $$$   if ~isempty(k)
           yy = y.(vn)(end-obs+1:end) + repmat(ys(i_var(j)),obs,1);
           plot(x(1:4),yy,'k-','LineWidth',1);
           hold on
% $$$       obs = 10;
% $$$   end
       plot(x,[NaN(obs-1,1); yy(end); yf.(vn)],'b-','LineWidth',1);
       hold on
       plot(x,[NaN(obs-1,1); yy(end); hpdinf.(vn)],'b--','LineWidth',1);
       hold on
       plot(x,[NaN(obs-1,1); yy(end); hpdsup.(vn)],'b--','LineWidth',1);
       hold on
       plot(x,[NaN(obs-1,1); yy(end); hpdinf2.(vn)],'b:','LineWidth',1);
       hold on
       plot(x,[NaN(obs-1,1); yy(end); hpdsup2.(vn)],'b:','LineWidth',1);
       title(vn,'Interpreter','none');
       hold off
       m = m + 1;
   end

   if m > 1
       eval(['print -deps ' fname '/graphs/forcst' int2str(n_fig) '.eps;'])
   end

