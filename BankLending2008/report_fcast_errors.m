% report_fcast_errors(fy, ndiffuse, fcast_steps, string_dates)
%
% makes graphs of n-step ahead fcast errors for given n's. The n's are
% given in a vector fcast_steps, fy, F, and ndiffuse should be a result of
% calc_fcast.
%
% maximum number in fcast_steps must be less or equal to size(fy,2)
%

function report_fcast_errors(fy, F, ndiffuse, fcast_steps, string_date)

global options_ M_

% set variables for which we want rmses of YoY and QoQ growths
grvars = strvcat('LGDP','LZ','LCPI','LGDP_EU','LCPI_EU','LCPI_US','LGDP_US', 'LGDP_EQ_US');

[first_year, first_quarter] = str_date_to_num(string_date);

eval(options_.datafile);

nvar = size(fy,1);
gend = size(fy,3)-1;

[xticks, xlabels] = yearly_ticks(first_year, first_quarter, gend, 8);

for h=fcast_steps
  disp(['----------- rmse ' num2str(h) ' ----------------------------------']);
  for ivar=1:nvar
    svar=deblank(options_.varobs(ivar,:));
    tit = ['Forecast error of ' svar ' ' int2str(h) ' periods ahead'];
    figure('Name', tit, 'Visible', 'off');
    title(tit);
    eval(['y=' svar ';']);
    y = y(1:gend);
    yfcast = reshape(fy(ivar,h,:), gend+1, 1);
    % do rmse
    rmse = sqrt(mean((y(ndiffuse+h:gend)-yfcast(ndiffuse+1:gend-h+1)).^2));
    ff = reshape(F(ivar,h,:), gend+1, 1);
    th_rmse = sqrt(mean(ff(ndiffuse+1:gend-h+1)));
    % printout
    disp(sprintf('%-15s %2d rmse=%10.6g th_rmse=%10.6g', svar, h, rmse, th_rmse));
    if ~isempty(strmatch(svar, grvars, 'exact'))
      % do QoQ
      if h == 1
        rmse_qoq = 4*rmse;
      else
        yfcastgr = 4*(reshape(fy(ivar,h,:),gend+1,1) - reshape(fy(ivar,h-1,:),gend+1,1));
        ygr = 4*(y(2:gend)-y(1:gend-1));
        rmse_qoq = sqrt(mean((ygr(ndiffuse+h-1:gend-1)-yfcastgr(ndiffuse+1:gend-h+1)).^2));
      end
      disp(sprintf('%-15s %2d rmse=%10.6g', [svar '(QoQ gr.)'], h, rmse_qoq));
      % do YoY
      if h <= 4
        rmse_yoy = sqrt(mean((y(ndiffuse+4:gend)-yfcast(ndiffuse+4-h+1:gend-h+1)).^2));
      else
        yfcastgr = reshape(fy(ivar,h,:),gend+1,1) - reshape(fy(ivar,h-4,:),gend+1,1);
        ygr = y(5:gend)-y(1:gend-4);
        rmse_yoy = sqrt(mean((ygr(ndiffuse+h-4+1:gend-4)-yfcastgr(ndiffuse+1:gend-h)).^2));
      end
      disp(sprintf('%-15s %2d rmse=%10.6g', [svar '(YoY gr.)'], h, rmse_yoy));      
    end
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

