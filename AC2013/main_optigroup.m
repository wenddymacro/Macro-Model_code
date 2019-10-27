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
countries   = {'Australia','Austria','Canada','Denmark','France','Iceland','Italy','Netherlands','Norway','Spain','Sweden','Switzerland','UK','USA'};
T0          = [1921 1947 1921 1835 1899 1838 1872 1850 1846 1908 1751 1876 1841 1933];
T0g         = [1920 1945 1920 1835 1895 1835 1870 1850 1845 1905 1750 1875 1840 1930];
T1          = [2004 2005 2004 2005 2004 2005 2003 2004 2005 2005 2005 2005 2003 2004];
age         = (0:110)';
nage        = size(age,1);
i60         = find(age==60);
i15         = find(age==15);
ngroup      = 4;
nbc         = size(T0,2);
scale       = [0.1 0.1 0.01 0.1 0.01 1 0.01 0.1 0.1 0.01 0.1 0.1 0.01 0.001];
for i=1:nbc;
    tic
    disp(countries{i})
    T   = T1(i)-T0(i)+1;
    date= (T0(i):T1(i))';
    file= ['data\pop_' char(countries{i}) ext];
    pop = load(file);
    
    popF= sum(reshape(pop(:,3),nage,T))';
    popM= sum(reshape(pop(:,4),nage,T))';
    popT= sum(reshape(pop(:,5),nage,T))';
    POP = [popT popM popF];
    pop = fliplr(pop(:,3:5));
    file= ['data\birth_' char(countries{i}) ext];
    brth= load(file);
    nat = fliplr(brth(:,2:4)./[popF popM popT]);
    file= ['data\death_' char(countries{i}) ext];
    dth = load(file);
    dth = fliplr(dth(:,3:5));
    P   = zeros(size(date,1),3);
    M   = zeros(size(date,1),3);
    A   = zeros(size(date,1),3);
    P0  = zeros(size(date,1),3);
    A0  = zeros(size(date,1),3);
    e60 = zeros(size(date,1),3);
    p60 = zeros(size(date,1),3);
    p15 = zeros(size(date,1),3);
    dr  = zeros(size(date,1),3);
    eva = zeros(size(date,1),3);
    amod= zeros(size(date,1),3);
    for j=1;%:3
        file= ['data\' char(countries{i}) '_' type{j} ext];
        data= load(file);
        eff = reshape(pop(:,j),nage,T);
        ev  = reshape(data(:,10),nage,T);
        tmp = reshape(dth(:,j),nage,T);
        mort= eff.*tmp;
        [m,im]=max(mort);
        amod(:,j)= age(im)';
        dr(:,j)= (sum(eff.*tmp)./sum(eff))';
        clear data;
        [p,a]=optimal_grouping_c(age,eff,ngroup,scale(i));
        p60(:,j) = (sum(eff(i60:nage,:))./sum(eff))';
        p15(:,j) = (sum(eff(1:i15,:))./sum(eff))';
        e60(:,j) = ev(i60,:)';
        P(:,j)=p(:,ngroup);
        A(:,j)=a(:,ngroup);
        P0(:,j)=p(:,1);
        A0(:,j)=a(:,2);
        M(:,j)=((age'*eff)./sum(eff))';
        for t=1:T;
            eva(t,j)=interp1(age,ev(:,t),A(t,j));
        end
    end
    eval(['save data\' char(countries{i}) ' date P0 A0 P A M p60 e60 eva nat dr pop amod;'])
end
