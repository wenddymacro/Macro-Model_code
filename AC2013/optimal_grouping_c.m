function [p,a]=optimal_grouping_c(age,eff,ngroup,scale)

[nobs,c]= size(eff);
p       = zeros(c,ngroup);
a       = zeros(c,ngroup+1);
agemin  = min(age);
agemax  = max(age);
for j=1:c;
    f   = eff(:,j)/sum(eff(:,j));
    tmp = [fix(eff(:,j)*scale);0];
    tmpage = [age:agemax+1];
    N   = sum(tmp);
    AGE = zeros(N,1);
    i0  = 1;    
    for i=1:(agemax-agemin+2);
        nb  = tmp(i);
        if nb>1
            AGE(i0:i0+nb-1)   = linspace(tmpage(i)+(tmpage(i+1)-tmpage(i))/(nb-1),tmpage(i+1),nb)';
            i0  = i0+nb;
        elseif nb==1
            AGE(i0)   = age(i);
            i0  = i0+1;
        end
    end
    a0      = linspace(agemin,agemax,ngroup+1);
    crit=1;
    tol=1e-6;
    while crit>tol;
        for i=2:ngroup;
            Id      = find((AGE>=a0(i-1))&(AGE<a0(i+1)));
            a(j,i)  = mean(AGE(Id));
        end
        a(j,1)          = a0(1);
        a(j,ngroup+1)   = a0(end);
        crit            = max(abs(a(j,:)-a0));
        a0              = a(j,:);
    end
    for i=1:ngroup
        Id      = ((AGE>=a(j,i))&(AGE<a(j,i+1)));
        p(j,i)  = sum(Id)/N;
    end
end
