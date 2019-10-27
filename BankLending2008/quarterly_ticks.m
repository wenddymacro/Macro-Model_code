% This generates numerical ticks and lables for the ticks. We suppose
% that the first item in a graph is associated with the given start_year
% and start_quarter and has x=1. We suposee that there are n items in total. The
% function generates ticks and labels with the given maximum number of
% yearly ticks. The ticks are at every quarter years.

function [ticks, labels] = quarterly_ticks(start_year, start_quarter, n, ...
                                        max_ticks)



% decide what ticks will get whole years
y_incr = ceil(n/4/max_ticks);
y_ticks = [mod(5-start_quarter,4)+1 : y_incr*4 : n];
% generate ticks and labels
labels = [];
ticks = 1:n;
yy = start_year;
qq = start_quarter-1;
for t = ticks;
  if isempty(find(y_ticks == t))
    if qq == 0 | qq == 2
      labels = strvcat(labels, ['Q' int2str(qq+1)]);
    else
      labels = strvcat(labels, ' ');
    end
  else
    yystr = int2str(yy);
    labels = strvcat(labels, [yystr(3:4) 'Q' int2str(qq+1)]);
  end
  qq = qq + 1;
  if (qq == 4)
    yy = yy + 1;
    qq = 0;
  end
end;
