function report_smoothed_rmse(fx, Fxall, ndiffuse, vars, fcast_steps, ...
                              string_date)

global options_ M_ oo_;

[first_year, first_quarter] = str_date_to_num(string_date);

nvar = size(vars,1);
gend = size(fx,3)-1;

[xticks, xlabels] = yearly_ticks(first_year, first_quarter, gend, 8);

for h=fcast_steps
  disp(['----------- rmse ' num2str(h) [' ----------------------------------']);
  for ivar=1:nvar
    svar = deblank(vars(ivars,:));
  end
end