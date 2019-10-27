function h=recplotband(X,Y,UP_Y,LW_Y)

col={'k','k','k','k'};
styl={'-','--','-','-'};
MS=[0.1 0.1 0.1 5];
h=plot(X,Y);
set(h,'linestyle','none')
% for i=1:length(h);
%     set(h(i),'linewidth',1,'color',char(col{i}),'linestyle',styl{i},'marker','o','markersize',MS(i));
% end
hpatch = patch([X;flipud(X)],[LW_Y;flipud(UP_Y)],0.7*ones(1,3)); 
set(hpatch,'linestyle','none')
hold on

h=plot(X,Y);
for i=1:length(h);
    set(h(i),'linewidth',1.2,'color',char(col{i}),'linestyle',styl{i},'marker','o','markersize',MS(i));
end
set(gca,'Layer','top')
grid
box on;
hold off;