% adapted for Dynare version 4. February 2008. Michel Juillard
%
function report_fcast_dynamic_4(fy, F, Udecomp, ndiffuse, freq, nahead, exo_ord, ...
                              string_date,var_list)

global options_ M_

[first_year, first_quarter] = str_date_to_num(string_date);

nvar = size(fy,1);
gend = size(fy,3)-1;

% decide about the number of batches, target batch length of 25
nb = ceil((gend-ndiffuse)/50);
batch_len = floor((gend-ndiffuse)/nb);

% cycle through batches and do them
for ib=1:nb
  q = first_quarter-1+ndiffuse+(ib-1)*batch_len;
  ffq = mod(q,4)+1;
  ffy = first_year+floor(q/4);
  if ib < nb
    len = batch_len;
  else
    len = gend-ndiffuse-(ib-1)*batch_len;
  end
  draw_dynamic_graph(fy, F, Udecomp, nahead, exo_ord, ffy, ffq, (ib-1)* ...
                     batch_len+ndiffuse, len, ib, var_list);
end



% function

function draw_dynamic_graph(fy, F, Udecomp, nahead, exo_ord, first_year, ...
                            first_quarter, offset, len, id,var_list)
global options_ M_ oo_

nvar = size(fy,1);

eval(options_.datafile);

% colormap
cmap = [
         0         0    0.5625
         0         0    0.6250
         0         0    0.6875
         0         0    0.7500
         0         0    0.8125
         0         0    0.8750
         0         0    0.9375
         0         0    1.0000
         0    0.0625    1.0000
         0    0.1250    1.0000
         0    0.1875    1.0000
         0    0.2500    1.0000
         0    0.3125    1.0000
         0    0.3750    1.0000
         0    0.4375    1.0000
         0    0.5000    1.0000
         0    0.5625    1.0000
         0    0.6250    1.0000
         0    0.6875    1.0000
         0    0.7500    1.0000
         0    0.8125    1.0000
         0    0.8750    1.0000
         0    0.9375    1.0000
         0    1.0000    1.0000
    0.0625    1.0000    0.9375
    0.1250    1.0000    0.8750
    0.1875    1.0000    0.8125
    0.2500    1.0000    0.7500
    0.3125    1.0000    0.6875
    0.3750    1.0000    0.6250
    0.4375    1.0000    0.5625
    0.5000    1.0000    0.5000
    0.5625    1.0000    0.4375
    0.6250    1.0000    0.3750
    0.6875    1.0000    0.3125
    0.7500    1.0000    0.2500
    0.8125    1.0000    0.1875
    0.8750    1.0000    0.1250
    0.9375    1.0000    0.0625
    1.0000    1.0000         0
    1.0000    0.9375         0
    1.0000    0.8750         0
    1.0000    0.8125         0
    1.0000    0.7500         0
    1.0000    0.6875         0
    1.0000    0.6250         0
    1.0000    0.5625         0
    1.0000    0.5000         0
    1.0000    0.4375         0
    1.0000    0.3750         0
    1.0000    0.3125         0
    1.0000    0.2500         0
    1.0000    0.1875         0
    1.0000    0.1250         0
    1.0000    0.0625         0
    1.0000         0         0
    0.9375         0         0
    0.8750         0         0
    0.8125         0         0
    0.7500         0         0
    0.6875         0         0
    0.6250         0         0
    0.5625         0         0
    0.5000         0         0
];

freq = 1;
[xticks, xlabels] = quarterly_ticks(first_year, first_quarter, len+nahead-1, 4);

order_var = oo_.dr.order_var;
exo_nbr = M_.exo_nbr;
labs = remove_prefix(M_.exo_names);
cml = size(cmap,1);
ind = floor(1:cml/(exo_nbr+1):cml);
cm = cmap(ind,:);

for ivar=1:nvar
    %    svar=options_.varobs(ivar,:);
    svar = M_.endo_names(order_var(ivar),:);
    if isempty(strmatch(svar,var_list,'exact'))
        continue
    end
  tit = [deblank(svar) ' ' num2str(nahead) ' periods ' ...
         'ahead (' num2str(id) ')'];
  figure('Name', tit, 'Visible', 'off');
%  eval(['y=' deblank(svar) ';']);
  y = oo_.SmoothedVariables.(deblank(svar))+oo_.dr.ys(order_var(ivar));
  % calculate th_rmse
  ff = F(ivar,nahead,:);
  ff = reshape(ff, size(ff,3), 1);
  th_rmse = sqrt(mean(ff(1:end-nahead)));

  % first plot with observable and forecasts
  subplot(3,1,1);
  rng = 1+offset:offset+len;
  plot(rng, y(rng), 'LineWidth', 3, 'Color', [0.7 0.7 0.7]);
  title(tit,'Interpreter','none');
  hold on
  for t = rng
    ycast = reshape(fy(ivar,1:nahead,t+1), nahead, 1);
    ycast = [y(t); ycast];
    plot(t:t+nahead, ycast);
    hold on
  end
  grid;
  set(gca,'XTickMode','manual');
  set(gca,'XTickLabelMode','manual');
  set(gca,'XTick',xticks+offset);
  set(gca,'XTickLabel',xlabels);
  set(gca,'XLimMode','manual');
  set(gca,'XLim', [offset offset+len+nahead]);
%   % second plot with forecast updates
%   subplot(3,1,2);
%   plot([1 len+nahead-1], [0 0]);
%   hold on
%   uall = reshape(Udecomp(ivar,nahead,rng,:), size(rng,2), exo_nbr);
%   u = sum(uall,2);
%   u = [zeros(nahead-1,1); u];
%   bar(u, 'stack','k');
%   if (nahead == 1)
%     hold on
%     plot([1 len],[-2*th_rmse -2*th_rmse],'Color',[0 0 0]);
%     hold on
%     plot([1 len],[2*th_rmse 2*th_rmse],'Color',[0 0 0]);
%   end
%   set(gca,'XlimMode','manual');
%   set(gca,'Xlim', [0 len+nahead]);
%   set(gca,'YLimMode','manual');
%   if (nahead == 1)
%     set(gca,'YLim',[1.2*min(min(u),-2*th_rmse) 1.2*max(max(u),+2* ...
%                                                       th_rmse)]);
%   else
%     set(gca,'YLim',[1.2*min(u) 1.2*max(u)]);
%   end
%   set(gca,'XTickMode','manual');
%   set(gca,'XTickLabelMode','manual');
%   set(gca,'XTick',xticks);
%   set(gca,'XTickLabel',xlabels);
% %   
%   % third plot with forecast updates decomposition
%   subplot(3,1,3);
%   [uall llabels sel] = select_large_u(uall, exo_ord, labs, 2);
%   ipos = find(uall>=0);
%   upos = zeros(size(uall));
%   upos(ipos) = uall(ipos);
%   ineg = find(uall<0);
%   uneg = zeros(size(uall));
%   uneg(ineg) = uall(ineg);
%   upos = [zeros(nahead-1,size(upos,2)); upos];
%   uneg = [zeros(nahead-1,size(uneg,2)); uneg];
%   bar(upos, 'stack');
%   hold on;
%   bar(uneg, 'stack');
%   colormap(cm([sel size(cm,1)],:));
%   set(gca,'XlimMode','manual');
%   set(gca,'Xlim', [0 len+nahead]);
%   set(gca,'YLimMode','manual');
%   set(gca,'YLim', [1.2*min(sum(uneg,2)) 1.2*max(sum(upos,2))]);
%   set(gca,'XTickMode','manual');
%   set(gca,'XTickLabelMode','manual');
%   set(gca,'XTick',xticks);
%   set(gca,'XTickLabel',xlabels);
%   hh = legend(llabels);
%   set(hh, 'Interpreter', 'none');
%   set(hh, 'Location', 'SouthOutside');
%   set(hh, 'Orientation','horizontal');
%   set(hh, 'FontSize',6);
%   set(gca,'XTickMode','manual');
%   set(gca,'XTickLabelMode','manual');
%   set(gca,'XTick',xticks);
%   set(gca,'XTickLabel',xlabels);
  % print
  print('-dpdf',[M_.dname '_dynamic_fcast_' deblank(svar) '_' num2str(nahead) '_' ...
                 num2str(id) '.pdf']);
  print('-depsc2',[M_.dname '_dynamic_fcast_' deblank(svar) '_' num2str(nahead) '_' ...
                   num2str(id) '.eps']);
end

function strout = remove_prefix(strin)

strout = [];
for i=1:size(strin,1)
  j = strfind(strin(i,:), '_');
  if isempty(j)
    strout = strvcat(strout, strin(i,:));
  else
    strout = strvcat(strout, strin(i,(j(1)+1):end));
  end
end


% this selects from given uall m-1 largest decompositions, the rest is
% summed to m-th, orders them according to exo_ord putting rest at the
% last place and returns the labels; sel is a vector of column indices of
% uallold corresponding to selected columns (its lenght is m-1)
function [uall, llabels sel] = select_large_u(uallold, exo_ord, labs, m)

if m > size(uallold,2)
  uall = uallold;
  llabels = labs(exo_ord,:);
  return
end

uall2 = uallold.^2;
uall2 = sum(uall2,1);
[tmp isort] = sort(uall2, 2, 'descend');

sel = [];
nonsel = [];
for iexo = exo_ord
  ff = find(isort == iexo);
  if ff >= m
    nonsel = [nonsel iexo];
  else
    sel = [sel iexo];
  end
end

uall = zeros(size(uallold,1),m);
uall(:,1:m-1) = uallold(:,sel);
uall(:,m) = sum(uallold(:,nonsel), 2);
llabels = labs(sel,:);
llabels = strvcat(llabels, 'REST');
