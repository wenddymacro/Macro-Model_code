clear
load shocktype
load resdisplay

if token == 1

figure;
subplot(4,2,1)
plot(oo_benchmark.irfs.gdp_eps_s*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.gdp_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.gdp_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.gdp_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Real GDP','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,2)
plot(oo_benchmark.irfs.pi_eps_s*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.pi_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.pi_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.pi_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Inflation','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,3)
plot(oo_benchmark.irfs.m2g_eps_s*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.m2g_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.m2g_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.m2g_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Money growth rate','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,4)
plot(((oo_benchmark.irfs.ltau_eps_s))*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(((oo_r.irfs.ltau_eps_s))*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(((oo_tau.irfs.ltau_eps_s))*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(((oo_both.irfs.ltau_eps_s))*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Required reserve ratio','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,5)
plot(oo_benchmark.irfs.ys_eps_s*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.ys_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.ys_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.ys_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('SOE output','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,6)
plot(oo_benchmark.irfs.yp_eps_s*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.yp_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.yp_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.yp_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('POE output','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,7)
plot(oo_benchmark.irfs.bcs_eps_s*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.bcs_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.bcs_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.bcs_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('SOE bankruptcy ratio','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,8)
plot(oo_benchmark.irfs.bcp_eps_s*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.bcp_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.bcp_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.bcp_eps_s*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('POE bankruptcy ratio','FontName','Times New Roman');
ylabel('percent deviations')
legend('Benchmark','Optimal money supply rule', 'Optimal \tau rule', 'Jointly optimal rules','Location','BestOutside');
suptitle('Impulse responses to TFP shock');
print -dpdf Figures/outfig_irfs_tfp_all.pdf;

end

if token == 2 
figure;
subplot(4,2,1)
plot(oo_benchmark.irfs.gdp_eps_as*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.gdp_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.gdp_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.gdp_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Real GDP','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,2)
plot(oo_benchmark.irfs.pi_eps_as*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.pi_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.pi_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.pi_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Inflation','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,3)
plot(oo_benchmark.irfs.m2g_eps_as*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.m2g_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.m2g_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.m2g_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Money growth rate','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,4)
plot(((oo_benchmark.irfs.ltau_eps_as))*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(((oo_r.irfs.ltau_eps_as))*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(((oo_tau.irfs.ltau_eps_as))*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(((oo_both.irfs.ltau_eps_as))*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Required reserve ratio','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,5)
plot(oo_benchmark.irfs.ys_eps_as*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.ys_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.ys_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.ys_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('SOE output','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,6)
plot(oo_benchmark.irfs.yp_eps_as*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.yp_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.yp_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.yp_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('POE output','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,7)
plot(oo_benchmark.irfs.bcs_eps_as*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.bcs_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.bcs_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.bcs_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('SOE bankruptcy ratio','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,8)
plot(oo_benchmark.irfs.bcp_eps_as*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.bcp_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.bcp_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.bcp_eps_as*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('POE bankruptcy ratio','FontName','Times New Roman');
ylabel('percent deviations')
legend('Benchmark','Optimal money supply rule', 'Optimal \tau rule', 'Jointly optimal rules','Location','BestOutside');
suptitle('Impulse responses to SOE TFP shock');
print -dpdf Figures/outfig_irfs_as_all.pdf;

end

if token == 3
figure;
subplot(4,2,1)
plot(oo_benchmark.irfs.gdp_eps_ap*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.gdp_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.gdp_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.gdp_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Real GDP','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,2)
plot(oo_benchmark.irfs.pi_eps_ap*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.pi_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.pi_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.pi_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Inflation','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,3)
plot(oo_benchmark.irfs.m2g_eps_ap*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.m2g_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.m2g_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.m2g_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Money growth rate','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,4)
plot(((oo_benchmark.irfs.ltau_eps_ap))*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(((oo_r.irfs.ltau_eps_ap))*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(((oo_tau.irfs.ltau_eps_ap))*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(((oo_both.irfs.ltau_eps_ap))*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('Required reserve ratio','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,5)
plot(oo_benchmark.irfs.ys_eps_ap*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.ys_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.ys_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.ys_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('SOE output','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,6)
plot(oo_benchmark.irfs.yp_eps_ap*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.yp_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.yp_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.yp_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('POE output','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,7)
plot(oo_benchmark.irfs.bcs_eps_ap*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.bcs_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.bcs_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.bcs_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('SOE bankruptcy ratio','FontName','Times New Roman');
ylabel('percent deviations')
subplot(4,2,8)
plot(oo_benchmark.irfs.bcp_eps_ap*100,'LineStyle','-','LineWidth',1.5,'Color',[0,0,0]);
hold all
plot(oo_r.irfs.bcp_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[0,0,1]);
hold all
plot(oo_tau.irfs.bcp_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,0]);
hold all
plot(oo_both.irfs.bcp_eps_ap*100,'LineStyle','--','LineWidth',1.5,'Color',[1,0,1]);
title('POE bankruptcy ratio','FontName','Times New Roman');
ylabel('percent deviations')
legend('Benchmark','Optimal money supply rule', 'Optimal \tau rule', 'Jointly optimal rules','Location','BestOutside');
suptitle('Impulse responses to POE TFP shock');
print -dpdf Figures/outfig_irfs_ap_all.pdf;

end
