% Pop_pays
%
%   Year          Age         Female    Male    Total
%
clear all
clc
ext         = '.txt';
type        = {'total','male','female'};
countries   = {'Australia','Austria','Canada','Denmark','France','Iceland','Italy','Netherlands','Norway','Spain','Sweden','Switzerland','UK','usa'};
T0          = [1921 1947 1921 1835 1899 1838 1872 1850 1846 1908 1751 1876 1841 1933];
T0g         = [1920 1945 1920 1835 1895 1835 1870 1850 1845 1905 1750 1875 1840 1930];
T1          = [2004 2005 2004 2005 2004 2005 2003 2004 2005 2005 2005 2005 2003 2004];
age         = (0:110)';
nage        = size(age,1);
nbc         = size(T0,2);
scale       = [0.001 0.001 0.001 0.01 0.001 0.1 0.001 0.01 0.01 0.001 0.01 0.01 0.001 0.001];

i           = 14;
disp(countries{i})
T           = T1(i)-T0(i)+1;
date        = (T0(i):T1(i))';
file        = ['data/pop_' char(countries{i}) ext];
pop         = load(file);
pop         = fliplr(pop(:,3:5));
eff         = reshape(pop(:,1),nage,T);
P           = zeros(T,5);
A           = zeros(T,5);
P0          = zeros(T,5);
A0          = zeros(T,5);
k           = 1;
for ngroup=2:6;
    [p,a]   = optimal_grouping_c(age,eff,ngroup,scale(i));
    P(:,k)  = p(:,ngroup);
    A(:,k)  = a(:,ngroup);
    P0(:,k)  = p(:,1);
    A0(:,k)  = a(:,2);
    k       = k+1;
end
Ecr = (P(:,2:end)./P0(:,2:end)); % Elder/Child ratio
% Generates Table 1
tmp = corrcoef(Ecr);
disp([(3:6)' tmp])
%%
fs=12;
smpl=1933:2004;

subplot(221);h=plot(smpl,P(:,2:end));
set(h(1),'linewidth',2,'color','k','linestyle','--');
set(h(2),'linewidth',2,'color','k','linestyle','-');
set(h(3),'linewidth',2,'color',0.5*ones(1,3),'linestyle','-');
set(h(4),'linewidth',2,'color',0.5*ones(1,3),'linestyle','--');
set(gca,'xlim',[1930 2005],'fontname','times','fontsize',fs);
xlabel('Years','fontname','times','fontsize',fs);
title('Share of Oldest Individuals','fontname','times','fontsize',fs);
grid

subplot(222);h=plot(smpl,P0(:,2:end));
set(h(1),'linewidth',2,'color','k','linestyle','--');
set(h(2),'linewidth',2,'color','k','linestyle','-');
set(h(3),'linewidth',2,'color',0.5*ones(1,3),'linestyle','-');
set(h(4),'linewidth',2,'color',0.5*ones(1,3),'linestyle','--');
xlabel('Years','fontname','times','fontsize',fs);
set(gca,'xlim',[1930 2005],'fontname','times','fontsize',fs);
title('Share of Youngest Individuals','fontname','times','fontsize',fs);
grid

subplot(212);
h=plot(smpl,Ecr);
set(h(1),'linewidth',2,'color','k','linestyle','--');
set(h(2),'linewidth',2,'color','k','linestyle','-');
set(h(3),'linewidth',2,'color',0.5*ones(1,3),'linestyle','-');
set(h(4),'linewidth',2,'color',0.5*ones(1,3),'linestyle','--');
xlabel('Years','fontname','times','fontsize',fs);
title('Elder-Child Ratio','fontname','times','fontsize',fs);
legend('n=3','n=4','n=5','n=6','location','southeast')
set(gca,'xlim',[1930 2005],'fontname','times','fontsize',fs);
grid

pause;
print('-dpdf','figure9');
close