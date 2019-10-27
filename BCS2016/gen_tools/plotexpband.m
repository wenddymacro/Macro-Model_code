function h=plotexpband(xt,xp,ssx,dp,tistr,varargin)

if nargin>5
    inver   = 0;
else
    inver   = 1;
end
T   = size(xp,1);

fs      = 12;
one     = ones(T,1);
xp0     = NaN*one;
xp1     = NaN*one;
id0     = find(dp<1);
id1     = find(dp>0.9);
if ~isempty(id0)
    xp0(id0)= xp(id0,1);
end
if ~isempty(id1)
    xp1(id1)= xp(id1,1);
end

hpatch = patch([xt;flipud(xt)],[xp(:,2);flipud(xp(:,3))],[.7 .7 .7]);
set(hpatch,'linestyle','none')
hold on

if inver
    if length(ssx)==1
        XG      = [xp(:,1) xp0 xp1 ssx*one];
    else
        XG      = [xp(:,1) xp0 xp1 ssx];
    end
    h       = plot(xt,XG);
    set(h(1),'linewidth',1,'color','k','linestyle',':');
    set(h(2),'linewidth',2,'color','k');
    set(h(3),'color','k','marker','o','markersize',3,'markerfacecolor','k');
    set(h(4),'linewidth',1,'color','k','linestyle','--');
    if length(h)==5
    set(h(5),'linewidth',1,'color','k','linestyle','--');
    end
else
    if length(ssx)==1
        XG      = [xp0 ssx*one];
    else
        XG      = [xp0 ssx];
    end
    h       = plot(xt,XG);
    set(h(1),'linewidth',1,'color','k');
    set(h(2),'linewidth',1,'color','k','linestyle','--');
end
set(gca,'fontname','times','fontsize',fs)
xlabel('Years','fontname','times','fontsize',fs)
title(tistr,'fontname','times','fontsize',fs)
box on;