function [year, quarter] = str_date_to_num(str_date);

%function start_d = str_date_to_num(str_date);
% Turns string dates to number dates.
% Input matrix must be in this form: 
%    str_date = ['1990Q1']
%
%Out put is in following form: 1998.0 - Quarter 1 1998
%                              2000.25 - Quarter 2 2000
%                              2002.5  - Quarter 3 2002
%                              2005.75 - Quarter 4 2005

year = str2num(str_date(1,1:4));
quarter = str2num(str_date(1,6));

