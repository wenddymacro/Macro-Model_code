fname = M_.dname;
% canada2;
% nobs = 51;

[xticks, xlabels] = yearly_ticks(1994, 1, nobs, 10);

%% observed variables report



figure;
plot(1:nobs,RS_US,1:nobs,oo_.SmoothedVariables.RR_US(1:nobs));
legend('RS\_US','RR\_US');
title('Observed Interest Rates');
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
print('-dpdf',[fname '_obs1.pdf']);
print('-depsc2',[fname '_obs1.eps']);


figure;
plot(5:nobs,oo_.SmoothedVariables.PIE_US(5:nobs)+pietar_us_ss, 5:nobs,oo_.SmoothedVariables.PIE_US4(5:nobs)+pietar_us_ss);
legend('PIE\_US', 'PIE\_US4');
title('Observed CPI\_US Inflation');
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
print('-dpdf',[fname '_obs3.pdf']);
print('-depsc2',[fname '_obs3.eps']);



%% equations' fit report
figure('Name','US Phillips curve');
title('US Phillips curve fit');
subplot(3,1,1);
plot(2:nobs, [oo_.SmoothedVariables.PIE_US(2:nobs)+pietar_us_ss ...
              oo_.SmoothedVariables.PIE_US(2:nobs)+pietar_us_ss-oo_.SmoothedShocks.RES_PIE_US(2:nobs)],...
              'LineWidth',2);
a1=gca;
legend('PIE\_US', 'PIE\_US Fit',...
       'Location','SouthOutside','Orientation','Horizontal');
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
subplot(3,1,2);
plot(1:nobs, oo_.SmoothedShocks.RES_PIE_US(1:nobs),'LineWidth',2);
a2=gca;
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
subplot(3,1,3);
a=zeros(nobs, 3);
a(2:nobs,1) = lambda_us1*(oo_.SmoothedVariables.E4_PIE_US4(2:nobs)+pietar_us_ss);
a(2:nobs,2) = (1-lambda_us1)*(oo_.SmoothedVariables.PIE_US4(1:nobs-1)+pietar_us_ss);
a(2:nobs,3) = lambda_us2*oo_.SmoothedVariables.Y_US(1:nobs-1);
ipos = find(a > 0);
apos = zeros(size(a));
apos(ipos) = a(ipos);
ineg = find(a < 0);
aneg = zeros(size(a));
aneg(ineg) = a(ineg);
bar(apos,'stack'); hold on; bar(aneg,'stack');
a3=gca;
hh=legend('US PIE4(+4) YoY', 'US PIE4(-1) YoY', 'Y\_US(-1)', ...
          'Location','SouthOutside','Orientation','Horizontal');
set(hh, 'FontSize',6);
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
%resize_subplot3(a1,a2,a3,0.40,0.20);
print('-dpdf',[fname '_PIEUS_fit.pdf']);
print('-depsc2',[fname '_PIEUS_fit.eps']);









figure('Name','US IS curve');
title('US IS curve fit');
subplot(3,1,1);
plot(1:nobs, [oo_.SmoothedVariables.Y_US(1:nobs) ...
              oo_.SmoothedVariables.Y_US(1:nobs) - oo_.SmoothedShocks.RES_Y_US(1:nobs)],...
     'LineWidth',2);
a1=gca;
legend('Y\_US', 'Y\_US Fit',...
       'Location','SouthOutside','Orientation','Horizontal');
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
subplot(3,1,2);
plot(1:nobs, oo_.SmoothedShocks.RES_Y_US(1:nobs),'LineWidth',2);
a2=gca;
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
subplot(3,1,3);
a=zeros(nobs, 3);
a(2:nobs,1) = beta_us1*oo_.SmoothedVariables.Y_US(1:nobs-1);
a(2:nobs,2) = beta_us2*oo_.SmoothedVariables.E1_Y_US(2:nobs);
a(2:nobs,3) = -beta_us3*(oo_.SmoothedVariables.RR_US(1:nobs-1)-oo_.SmoothedVariables.RR_US_BAR(1:nobs-1));
a(2:nobs,4) = -theta*oo_.SmoothedVariables.E2(2:nobs);
ipos = find(a > 0);
apos = zeros(size(a));
apos(ipos) = a(ipos);
ineg = find(a < 0);
aneg = zeros(size(a));
aneg(ineg) = a(ineg);
bar(apos,'stack'); hold on; bar(aneg,'stack');
a3=gca;
hh=legend('Y\_US(-1)', 'Y\_US(+1)', 'RR\_US\_GAP(-1)', 'RES\_BLT\_US',...
          'Location','SouthOutside','Orientation','Horizontal');
set(hh, 'FontSize',6);
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
%resize_subplot3(a1,a2,a3,0.40,0.20);
print('-dpdf',[fname '_YUS_fit.pdf']);
print('-depsc2',[fname '_YUS_fit.eps']);


figure('Name','US Interest rate rule');
title('US RS rule fit');
subplot(3,1,1);
plot(1:nobs, [oo_.SmoothedVariables.RS_US(1:nobs)+rr_us_bar_ss+2 ...
              oo_.SmoothedVariables.RS_US(1:nobs)+rr_us_bar_ss+2-oo_.SmoothedShocks.RES_RS_US(1:nobs)],...
              'LineWidth',2);
a1=gca;
legend('RS\_US', 'RS\_US Fit',...
       'Location','SouthOutside','Orientation','Horizontal');
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
subplot(3,1,2);
plot(1:nobs, oo_.SmoothedShocks.RES_RS_US(1:nobs),'LineWidth',2);
a2=gca;
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
subplot(3,1,3);
a=zeros(nobs, 4);

a(2:nobs,1) = gamma_us1*(oo_.SmoothedVariables.RS_US(1:nobs-1)+rr_us_bar_ss+2);
a(2:nobs,2) = (1-gamma_us1)*(oo_.SmoothedVariables.RR_US_BAR(2:nobs)+rr_us_bar_ss+2+oo_.SmoothedVariables.E4_PIE_US4(2:nobs));
a(2:nobs,3) = (1-gamma_us1)*gamma_us2*(oo_.SmoothedVariables.E4_PIE_US4(2:nobs)-pietar_us_ss);
a(2:nobs,4) = (1-gamma_us1)*gamma_us4*(oo_.SmoothedVariables.Y_US(2:nobs));

ipos = find(a > 0);
apos = zeros(size(a));
apos(ipos) = a(ipos);
ineg = find(a < 0);
aneg = zeros(size(a));
aneg(ineg) = a(ineg);
bar(apos,'stack'); hold on; bar(aneg,'stack');
a3=gca;
hh=legend('RS\_US(-1)', 'RS\_US neutral', 'Dev from pietar', 'Y\_US',...
          'Location','SouthOutside','Orientation','Horizontal');
set(hh, 'FontSize',6);
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
%resize_subplot3(a1,a2,a3,0.40,0.20);
print('-dpdf',[fname '_RSUS_fit.pdf']);
print('-depsc2',[fname '_RSUS_fit.eps']);

figure('Name','GDP_US');
title('GDP\_US');
plot(1:nobs, [oo_.SmoothedVariables.LGDP_US(1:nobs)+(1:nobs)'*growth_us_ss/4 ...
              oo_.SmoothedVariables.LGDP_US_BAR(1:nobs)+(1:nobs)'*growth_us_ss/4]);
legend('Obs. LGDP\_US','LGDP\_US\_BAR','Location','SouthOutside','Orientation','Horizontal');
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
print('-dpdf',[fname '_LGDPUS.pdf']);
print('-depsc2',[fname '_LGDPUS.eps']);

figure('Name','DOT_LGDP_US');
title('DOT\_LGDP\_US');
plot(2:nobs, [4*(LGDP_US(2:nobs)-LGDP_US(1:nobs-1)) ...
              4*(oo_.SmoothedVariables.LGDP_US_BAR(2:nobs)- ...
                   oo_.SmoothedVariables.LGDP_US_BAR(1:nobs-1))+growth_us_ss]);
legend('Obs. DOT\_LGDP\_US','DOT\_LGDP\_US\_BAR','Location','SouthOutside','Orientation','Horizontal');
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
print('-dpdf',[fname '_LGDPUS_BAR.pdf']);
print('-depsc2',[fname '_LGDPUS_BAR.eps']);


figure('Name','UNR_US');
title('UNR\_US');
plot(1:nobs, [UNR_US UNR_US+oo_.SmoothedVariables.UNR_US_GAP(1:nobs)]);
legend('Obs. UNR\_US','UNR\_US\_BAR','Location','SouthOutside','Orientation','Horizontal');
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
print('-dpdf',[fname '_UNRUS.pdf']);
print('-depsc2',[fname '_UNRUS.eps']);

figure('Name', 'UNR_US_GAP');
title('UNR\_US\_GAP');
plot(1:nobs, [oo_.SmoothedVariables.UNR_US_GAP(1:nobs) ...
              oo_.SmoothedVariables.Y_US(1:nobs)]);
legend('UNR\_US\_GAP', 'Y\_US','Location','SouthOutside','Orientation','Horizontal');
set(gca,'XTickMode','manual');
set(gca,'XTickLabelMode','manual');
set(gca,'XTick',xticks);
set(gca,'XTickLabel',xlabels);
grid;
print('-dpdf', [fname '_UNRUS_GAP.pdf']);
print('-depsc2', [fname '_UNRUS_GAP.eps']);

