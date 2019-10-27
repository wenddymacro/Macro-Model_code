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
file        = 'data/pop_usa.txt';
pop         = load(file);
T0          = 1933;
T1          = 2004;
age         = (0:110)';
nage        = size(age,1);
i60         = find(age==60);
i15         = find(age==15);
i65         = find(age==65);
gray        = 0.5*ones(1,3);
T           = T1-T0+1;
date        = (T0:T1)';

POP         = reshape(pop(:,5),nage,T);
POPA        = repmat(age,1,T).*POP;
CPA         = cumsum(POPA)./repmat(sum(POPA),111,1);
CP          = cumsum(POP)./repmat(sum(POP),111,1);
FPOP        = POP./repmat(sum(POP),111,1);

%%
t0          = find(date==1950);
t1          = find(date==2000);
fs          = 12;
h=stairs(100*FPOP(:,[t0 t1]),[age age]);
set(h(1),'color',0.65*ones(1,3),'linewidth',2);
set(h(2),'color','k','linewidth',2);
legend('1950','2000','location','Northeast');
xlabel('% of Total Population','fontname','times','fontsize',16)
ylabel('Age','fontname','times','fontsize',16)
set(gca,'ylim',[0 110],'fontsize',16,'fontname','times');
grid
print('-dpdf','figure2')
%%
h=plot(CP(:,[t0 t1]),CPA(:,[t0 t1]));
set(h(1),'color','k','linewidth',2,'linestyle','--');
set(h(2),'color','k','linewidth',2);
legend('1950','2000','location','NorthWest')
line([0 1],[0 1],'color','k')
xlabel('Cdf of Total Population','fontname','times','fontsize',16)
ylabel('Cdf of Total Lived Years','fontname','times','fontsize',16)
set(gca,'fontsize',16,'fontname','times');
print('-dpdf','figure3')
