function h=plotirf(xp,sx,varargin)

T       = length(sx);
hpatch  = patch([(1:T)';flipud((1:T)')],[xp(:,1)-sx;flipud(xp(:,1)+sx)],[.7 .7 .7]); 
set(hpatch,'linestyle','none')
hold on
one     = ones(T,1);
fs      = 12;
if nargin>2
    XG      = [xp varargin{1}*one];
    h       = plot(XG);
    set(h(1),'linewidth',1.5,'color','k','linestyle','-');
    set(h(2),'linewidth',1.5,'color','k','linestyle','--');
    set(h(3),'linewidth',1,'color','k','linestyle','--');
else
    h       = plot(xp);
    set(h(1),'linewidth',1.5,'color','k','linestyle','-');
    set(h(2),'linewidth',1.5,'color','k','linestyle','--');
end
set(gca,'fontname','times','fontsize',fs)
xlabel('Years','fontname','times','fontsize',fs)
box on;