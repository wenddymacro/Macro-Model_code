% This generates numerical ticks and lables for the ticks. We suppose
% that the first item in a graph is associated with the given start_year
% and start_quarter and has x=1. We suposee that there are n items in total. The
% function generates ticks and labels with the given maximum number of
% ticks. The ticks are at whole years.

function [ticks, labels] = yearly_ticks(start_year, start_quarter, n, ...
                                        max_ticks)



year = ceil(start_year+(start_quarter-1)/4);
y_incr = ceil(n/4/max_ticks);
ticks = [mod(5-start_quarter,4)+1 : y_incr*4 : n];
labels = [];
y = year;
for t = ticks;
  labels = strvcat(labels, int2str(y));
  y = y + y_incr;
end;
