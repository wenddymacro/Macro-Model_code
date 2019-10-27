% adapted for Dynare version 4. February 2008. Michel Juillard
%
function report_smoothed_errors_4a(fx, Fxall, ndiffuse, vars, fcast_steps, ...
                                string_date)

global options_ M_ oo_

[first_year, first_quarter] = str_date_to_num(string_date);

nvar = size(vars,1);
gend = size(fx,3)-1;
lgy_ = M_.endo_names;
ys = oo_.dr.ys;
order_var = oo_.dr.order_var;
[xticks, xlabels] = yearly_ticks(first_year, first_quarter, gend, 8);

for h=fcast_steps
  %disp(['& RMSE ' num2str(h) 'Q ahead & Theoretical RMSE ' num2str(h) 'Q ahead\\']);
  disp(['& RMSE ' num2str(h) 'Q ahead\\ ']);
  disp(['& \hline \\']);
  for ivar=1:nvar
    svar = deblank(vars(ivar,:));
    i = strmatch(svar, lgy_(order_var,:),'exact');
    is = strmatch(svar, lgy_,'exact');
    tit = ['Smoothed forecast error of ' svar ' ' int2str(h) ' periods ahead'];
    figure('Name', tit, 'Visible', 'off');
    title(tit);
    y = oo_.SmoothedVariables.(svar)(1:gend)+ys(is);
    yfcast = reshape(fx(i,h,:),gend+1,1);
    % do rmse
    rmse = sqrt(mean((y(ndiffuse+h:gend)-yfcast(ndiffuse+1:gend-h+1)).^2));
        rep(h,ivar) = rmse;
        
    ff = reshape(Fxall(ivar,h,:), gend+1, 1);
    th_rmse = sqrt(mean(ff(ndiffuse+1:gend-h+1)));
    % printout
    % printout
    disp([sprintf('$%-15s$ & %10.6g', svar, rmse) '\\']);
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

  disp(['& 1 Q Ahead & 4 Q Ahead & 8 Q Ahead & 12 Q Ahead \\ ']);
  disp(['& \hline \\']);
for ivar=1:nvar    
    svar = deblank(vars(ivar,:));
    disp([sprintf('$%-15s$ & %10.2g & %10.2g & %10.2g & %10.2g', svar, rep(1,ivar), rep(4,ivar), rep(8,ivar), rep(12,ivar)) '\\']);
    
end

    
    
    
    
