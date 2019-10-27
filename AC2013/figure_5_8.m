% Pop_pays & death_pays
%
%   Year          Age         Female    Male    Total
%
% birth_pays
%
%   Year          Female    Male    Total
%
%
%   Year          Age         mx       qx    ax      lx      dx      Lx       Tx     ex
%
%   1   Year Year or range of years (for both period & cohort data)
%   2   Age  Age group for n-year interval from exact age x to just before exact age x+n, where n=1, 4, 5, or ? (open age interval)
%   3   m(x)
%   4   q(x) Probability of death between exact ages x and x+n
%   5   a(x)
%   6   l(x) Number of survivors at exact age x, assuming l(0) = 100,000
%   7   d(x) Number of deaths between exact ages x and x+n
%   8   L(x) Number of person-years lived between exact ages x and x+n
%   9   T(x) Number of person-years remaining after exact age x
%  10   e(x) Life expectancy at exact age x (in years)

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
i60         = find(age==60);
i15         = find(age==15);
i65         = find(age==65);
ngroup      = 4;
nbc         = size(T0,2);
gris        = 0.5*ones(1,3);
fs          = 12;
scale       = [0.1 0.1 0.01 0.1 0.01 1 0.01 0.1 0.1 0.01 0.1 0.1 0.01 0.001];
i           = 14;   % set to 1:nbc to compute all countries

disp(countries{i})
T           = T1(i)-T0(i)+1;
date        = (T0(i):T1(i))';
file        = ['data/pop_' char(countries{i}) ext];
pop         = load(file);
POP         = reshape(pop(:,5),nage,T);
FPOP        = POP./repmat(sum(POP),111,1);
pop         = fliplr(pop(:,3:5));
eff         = reshape(pop(:,1),nage,T);
[p,a]       = optimal_grouping_c(age,eff,ngroup,scale(i));

%%
T           = 1933:2004;
TT          = [find(T==1950);find(T==2000)];
for k=1:2
    t0=TT(k);
    subplot(2,2,k)
    for i=1:4;
        patch([0 100*p(t0,i)/(a(t0,i+1)-a(t0,i)) 100*p(t0,i)/(a(t0,i+1)-a(t0,i)) 0],[a(t0,i) a(t0,i) a(t0,i+1) a(t0,i+1)],0.85*ones(1,3))
    end
    hold on
    h=stairs(100*FPOP(:,t0),age);
    set(h(1),'color','k','linewidth',2);
    title(T(TT(k)),'fontname','times','fontsize',14)
    xlabel('% of Total Population','fontname','times','fontsize',14)
    ylabel('Age','fontname','times','fontsize',14)
    set(gca,'ylim',[0 110],'xlim',[0 2.5])
    set(gca,'fontsize',16,'fontname','times');
    line([0 2.5],[a(TT(1),4) a(TT(1),4)],'linestyle','--','color','k')
    box on
    hold on
    grid
end
print('-dpdf','figure5')
%%
TT          = [find(T==1933);find(T==2004)];
smpl        = TT(1):TT(2);
subplot(221);h=plot(T,a(smpl,4));
set(h,'linewidth',2,'color','k','linestyle','-');
set(gca,'xlim',[1930 2005],'fontname','times','fontsize',fs);
xlabel('Years','fontname','times','fontsize',fs);
ylabel('Age in Years','fontname','times','fontsize',fs);
title('3^{rd} Cutoff Age','fontname','times','fontsize',fs);
grid
subplot(222);h=plot(T,p(smpl,4));
set(h,'linewidth',2,'color','k','linestyle','-');
set(gca,'ylim',[0.175 0.225],'xlim',[1930 2005],'fontname','times','fontsize',fs);
xlabel('Years','fontname','times','fontsize',fs);
title('Share of the Oldest Individuals','fontname','times','fontsize',fs);
grid
pause;
print('-dpdf','figure6');
close
%%
subplot(221);
h=plot(T,a(smpl,2));
set(h,'linewidth',2,'color','k','linestyle','-');
set(gca,'xlim',[1930 2005],'fontname','times','fontsize',fs);
xlabel('Years','fontname','times','fontsize',fs);
ylabel('Age in Years','fontname','times','fontsize',fs);
title('1^{st} Cutoff Age','fontname','times','fontsize',fs);
grid
subplot(222);
h=plot(T,p(smpl,1));
set(h,'linewidth',2,'color','k','linestyle','-');
set(gca,'ylim',[0.225 0.3],'xlim',[1930 2005],'fontname','times','fontsize',fs);
xlabel('Years','fontname','times','fontsize',fs);
title('Share of the Youngest Individuals','fontname','times','fontsize',fs);
grid
pause;
print('-dpdf','figure7');
close
%%
h=plot(T,p(smpl,4)./p(smpl,1));
set(h,'linewidth',2,'color','k','linestyle','-');
set(gca,'xlim',[1930 2005],'fontname','times','fontsize',16);
xlabel('Years','fontname','times','fontsize',16);
grid
pause;
print('-dpdf','figure8');
close