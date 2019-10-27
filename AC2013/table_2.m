clear all
clc
ext         = '.txt';
type        = {'total','male','female'};
countries   = {'Australia','Austria','Canada','Denmark','France','Iceland','Italy','Netherlands','Norway','Spain','Sweden','Switzerland','UK','USA'};
countriesl  = {'Australia','Austria','Canada','Denmark','France','Iceland','Italy','Netherlands','Norway','Spain','Sweden','Switzerland','England & Wales','USA'};
T0          = [1921 1947 1921 1835 1899 1838 1872 1850 1846 1908 1751 1876 1841 1933];
T0g         = [1920 1945 1920 1835 1895 1835 1870 1850 1845 1905 1750 1875 1840 1930];
T1          = [2004 2005 2004 2005 2004 2005 2003 2004 2005 2005 2005 2005 2003 2004];
age         = (0:110)';
nage        = size(age,1);
i60         = find(age==60);
ngroup      = 4;
nbc         = size(T0,2);
gris        = 0.5*ones(1,3);
fs          = 12;
sel         = [1:7 9 11:14];
S=[];
STD=[];
for ii=1:12;
    i       = sel(ii);
    T       = T1(i)-T0(i)+1;
    nm      = floor(0.05*T);
    date    = (T0(i):T1(i))';
    i0      = find(date==T1(i)-50);
    i1      = find(date==T1(i));
    file    = ['data/' char(countries{i})];
	load(file);
    STD=[STD std(P0(:,1))];
    reso=ols(log(P(:,1)),[ones(T,1) (1:T)']);
    reso50=ols(log(P(i0:i1,1)),[ones(i1-i0+1,1) (1:i1-i0+1)']);
    resy=ols(log(P0(:,1)),[ones(T,1) (1:T)']);
    resy50=ols(log(P0(i0:i1,1)),[ones(i1-i0+1,1) (1:i1-i0+1)']);
    resr=ols(log(P(:,1)./P0(:,1)),[ones(T,1) (1:T)']);
    resr50=ols(log(P(i0:i1,1)./P0(i0:i1,1)),[ones(i1-i0+1,1) (1:i1-i0+1)']);
    S=strvcat(S,sprintf('%s & %6.4f & %6.4f & %6.4f \\\\',countries{i},100*reso50.beta(2),100*resy50.beta(2),100*resr50.beta(2)));
    S=strvcat(S,sprintf(' & {\\scriptsize [%6.4f]} & {\\scriptsize [%6.4f]} & {\\scriptsize [%6.4f]} \\\\', ...
    tdis_prb(reso50.tstat(2),reso50.nobs-reso50.nvar),tdis_prb(resy50.tstat(2),resy50.nobs-resy50.nvar),tdis_prb(resr50.tstat(2),resr50.nobs-resr50.nvar)));
end
disp(S)