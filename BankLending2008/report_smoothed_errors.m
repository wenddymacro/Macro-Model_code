function report_smoothed_errors(fx, Fxall, ndiffuse, vars, fcast_steps, ...
                                string_date)

global options_ M_ oo_ lgy_ dr_;

[first_year, first_quarter] = str_date_to_num(string_date);

nvar = size(vars,1);
gend = size(fx,3)-1;

[xticks, xlabels] = yearly_ticks(first_year, first_quarter, gend, 8);

for h=fcast_steps
  disp(['----------- rmse ' num2str(h) ' ----------------------------------']);
  for ivar=1:nvar
    svar = deblank(vars(ivar,:));
    i = strmatch(svar, lgy_(dr_.order_var,:),'exact');
    tit = ['Smoothed forecast error of ' svar ' ' int2str(h) ' periods ahead'];
    figure('Name', tit, 'Visible', 'off');
    title(tit);
    y = oo_.SmoothedVariables.(svar)(1:gend);
    yfcast = reshape(fx(i,h,:),gend+1,1);
    % do rmse
    rmse = sqrt(mean((y(ndiffuse+h:gend)-yfcast(ndiffuse+1:gend-h+1)).^2));
    ff = reshape(Fxall(ivar,h,:), gend+1, 1);
    th_rmse = sqrt(mean(ff(ndiffuse+1:gend-h+1)));
    % printout
    disp(sprintf('%-15s %2d rmse=%10.6g th_rmse=%10.6g', svar, h, rmse, th_rmse));
    % first plot
    subplot(2, 1, 1);
    plot(1:gend, y, ndiffuse+h:gend, yfcast(ndiffuse+1:gend-h+1));
    set(gca,'XTickMode','manual');
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTick',xticks);
    set(gca,'XTickLabel',xlabels);
    grid;
    hh = legend(deblank(svar), [deblank(svar) '_fcast']);
    set(hh, 'Interpreter', 'none');
    set(hh, 'Location','NorthWest');
    % second plot
    subplot(2, 1, 2);
    plot([1 gend],[0 0],'-r');
    hold on
    plot(ndiffuse+h:gend, yfcast(ndiffuse+1:gend-h+1)-y(ndiffuse+h:gend),'-k');
    hold on
    plot([ndiffuse+h gend], [rmse rmse], '--');
    hold on
    plot([ndiffuse+h gend], [-rmse -rmse], '--');
    set(gca,'XTickMode','manual');
    set(gca,'XTickLabelMode','manual');
    set(gca,'XTick',xticks);
    set(gca,'XTickLabel',xlabels);
    grid;   
    hold off
    % print
    print('-dpdf',[M_.dname '_fcast_error_' deblank(svar) '_' num2str(h) ...
                   '.pdf']); 
    print('-depsc2',[M_.dname '_fcast_error_' deblank(svar) '_' num2str(h) ...
                     '.eps']); 
  end
end