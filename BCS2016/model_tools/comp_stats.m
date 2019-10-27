function comp_stats(file,table,gr_)
%%
load(strcat('results/',file,'_approx'));   % loads the file containing approximation data
load(strcat('results/',file,'_simul'));    % loads the file containing simulated data
psi     = 1.2;  % Growth rate
lbhp    = 6.25; % HP filter parameter (annual, see Uhlig and Ravn)
T       = length(sim.y);
dy      = [0;diff(log(sim.y))+psi/100];
hpy     = log(sim.y)-hpfilter(log(sim.y),lbhp);
hpk     = log(sim.k)-hpfilter(log(sim.k),lbhp);
xy      = log(sim.y)+psi*(1:T)'/100;
xk      = log(sim.k)+psi*(1:T)'/100;
xh      = log(sim.h);
xz      = log(sim.z);
%% Identify Banking Crises

ptc     = find(diff(sim.dum)==1);  % Find initial date before crisis (+1 to get the crisis itself)
etc     = find(diff(sim.dum)==-1); % Find last date of crisis
Nc      = length(ptc);
%% Identifies Recessions
% compute the dates of the recessions NOT associated with banking crises
% within the recessions (be it at the start or at the end of the crisis)
pc      = 27.2;             % set such that we get same % of recessions as in the data
dy0     = prctile(dy,pc)-psi/100;
Id      = find(dy<dy0);    % y(t)>y(t+1) => potential start of recession
ni      = length(Id);
dumr    = zeros(T,1);
dumr(Id)= 1;
for i=1:ni
    j   = Id(i);
    crit= 1;
    k   = 1;
    while and((crit==1),(j+k<=T-1))
        crit        = dy(j+k)<0;
        dumr(j+k)   = crit;
        k           = k+1;
    end
end
ddr     = [diff(dumr);0];          % 1=> recession starts next period, -1 => ends this period
rpeak   = find(ddr==1);            % Peak before recession
rtrough = find(ddr==-1);           % trough of recession
Nr      = min(length(rpeak),length(rtrough));
rpeak   = rpeak(1:Nr);
rtrough = rtrough(1:Nr);
% only keep the recessions with more than per quarter of negative growth
% as in the NBER, the recessions starts at the first period of negative
% growth (provided that at least one negative growth period is to come)
% (NBER recessions go from peaks to troughs)
per     = 0;
id      = find((rtrough-rpeak)>=per);
TT      = [rpeak(id) rtrough(id)];
Nr      = size(TT,1);

durr    = zeros(Nr,1);
dyptr   = zeros(Nr,1);
dhpyptr = zeros(Nr,1);
dzptr   = zeros(Nr,1);
dhptr   = zeros(Nr,1);
dy2r    = zeros(Nr,1);
dhpkptr = zeros(Nr,1);
dk2r    = zeros(Nr,1);
dkm2r   = zeros(Nr,1);
dhpk2r  = zeros(Nr,1);
dhpkm2r = zeros(Nr,1);
devk2r  = zeros(Nr,1);
hpyr    = zeros(Nc,13);
hpkr    = zeros(Nc,13);

for i=1:Nr;
    durr(i)     = rtrough(i)-rpeak(i);
    dyptr(i)    = 100*(xy(rtrough(i))-xy(rpeak(i)));
    dhpyptr(i)  = 100*(hpy(rtrough(i))-hpy(rpeak(i)));
    dzptr(i)    = 100*(xz(rtrough(i))-xz(rpeak(i)));
    dhptr(i)    = 100*(xh(rtrough(i))-xh(rpeak(i)));
    dy2r(i)     = 100*(xy(rpeak(i)+2)-xy(rpeak(i)));
    dhpkptr(i)  = 100*(hpk(rtrough(i))-hpk(rpeak(i)));
    dk2r(i)     = 100*(xk(rpeak(i)+2)-xk(rpeak(i)));
    dhpk2r(i)   = 100*(hpk(rpeak(i)+2)-hpk(rpeak(i)));
    dhpkm2r(i)  = 100*(hpk(rpeak(i))-hpk(rpeak(i)-2));
    devk2r(i)   = 100*hpk(rpeak(i));
    hpyr(i,:)   = 100*hpy(rpeak(i)-5:rpeak(i)+7);
    hpkr(i,:)   = 100*hpk(rpeak(i)-5:rpeak(i)+7);
end

%%
% we get all the recessions dates, including those that involve financial crises
recwobc = [];
dwo     = zeros(length(TT),1);
for i=1:length(TT);
    if ~any(ptc>=TT(i,1) & ptc<TT(i,2))
        recwobc = [recwobc;TT(i,:)];
        dwo(i)  = 1;
    end
end
Nrwo    = sum(dwo);
Nrw     = Nr-sum(dwo);
%%
dys     = prctile(dyptr,100/3);
dyl     = prctile(dyptr,200/3);
MDUR    = [mean(durr(dwo==0)) mean(durr(dwo==1)) mean(durr) mean(durr(dyptr<=dys)) mean(durr(dyptr>=dyl))];
MDHPYPT = [mean(dhpyptr(dwo==0)) mean(dhpyptr(dwo==1)) mean(dhpyptr) mean(dhpyptr(dyptr<=dys)) mean(dhpyptr(dyptr>=dyl))];
MDYPT   = [mean(dyptr(dwo==0)) mean(dyptr(dwo==1)) mean(dyptr) mean(dyptr(dyptr<=dys)) mean(dyptr(dyptr>=dyl))];
MDY2    = [mean(dy2r(dwo==0)) mean(dy2r(dwo==1)) mean(dy2r) mean(dy2r(dyptr<=dys)) mean(dy2r(dyptr>=dyl))];
MDHPKPT = [mean(dhpkptr(dwo==0)) mean(dhpkptr(dwo==1)) mean(dhpkptr) mean(dhpkptr(dyptr<=dys)) mean(dhpkptr(dyptr>=dyl))];
MDK2    = [mean(dk2r(dwo==0)) mean(dk2r(dwo==1)) mean(dk2r) mean(dk2r(dyptr<=dys)) mean(dk2r(dyptr>=dyl))];
MDKM2   = [mean(dkm2r(dwo==0)) mean(dkm2r(dwo==1)) mean(dkm2r) mean(dkm2r(dyptr<=dys)) mean(dkm2r(dyptr>=dyl))];
MDHPK2  = [mean(dhpk2r(dwo==0)) mean(dhpk2r(dwo==1)) mean(dhpk2r) mean(dhpk2r(dyptr<=dys)) mean(dhpk2r(dyptr>=dyl))];
MDHPKM2 = [mean(dhpkm2r(dwo==0)) mean(dhpkm2r(dwo==1)) mean(dhpkm2r) mean(dhpkm2r(dyptr<=dys)) mean(dhpkm2r(dyptr>=dyl))];
MDEVK2  = [mean(devk2r(dwo==0)) mean(devk2r(dwo==1)) mean(devk2r) mean(devk2r(dyptr<=dys)) mean(devk2r(dyptr>=dyl))];
%%
if table==1
    disp('\begin{table}[htbp]')
    disp('\centering')
    disp('\caption{Model Statistics on Recessions}\label{tab:modtabrec}')
    disp('\vspace{0.5\baselineskip}')
    disp('\begin{tabular}{lccccc}')
    disp('\toprule')
    disp('& Financial & Other &~& All &~& Severe & Mild \\')
    disp('\midrule')
    fprintf('Frequency (\\%%) & %4.2f & %4.2f && %4.2f && %4.2f & %4.2f \\\\ \n',100*[Nrw/T Nrwo/T Nr/T Nr/(3*T) Nr/(3*T)])
    fprintf('Duration (peak to trough, years)      & %6.2f & %6.2f && %6.2f && %6.2f & %6.2f \\\\ \n',MDUR)
    fprintf('Output Loss (Peak to trough, \\%%)      & %6.2f & %6.2f && %6.2f && %6.2f & %6.2f \\\\ \n',MDYPT)
    disp('\midrule')
    disp('\multicolumn{6}{l}{\textit{Recession (\% unless otherwise)}}\\')
    fprintf('HP--Credit Growth (Peak to trough, \\%%)    & %6.2f & %6.2f && %6.2f && %6.2f & %6.2f \\\\ \n',MDHPKPT)
    fprintf('HP--Credit Growth (First 2 years, \\%%) & %6.2f & %6.2f && %6.2f && %6.2f & %6.2f \\\\ \n',MDHPK2)
    disp('\midrule')
    disp('\multicolumn{6}{l}{\textit{Pre-Recession (\%)}}\\')
    fprintf('HP--Credit Growth (2 years before peak) & %6.2f & %6.2f && %6.2f && %6.2f & %6.2f \\\\ \n',MDHPKM2)
    fprintf('HP dev. of Credit at peak & %6.2f & %6.2f && %6.2f && %6.2f & %6.2f \\\\ \n',MDEVK2)
    disp('\bottomrule')
    disp('\end{tabular}')
    
    disp('\parbox{\textwidth}{\footnotesize\underline{Note:} }')
    disp('\end{table}')
else
    disp(Nr/T)
    fprintf('& %6.2f \n',mean(sim.rho(sim.dum==0)))
    fprintf('& %6.2f \n',mean(sim.R))
    fprintf('& %6.2f \n',mean(sim.r))
    fprintf('& %6.2f \n',mean(sim.R-sim.r))
    fprintf('& %6.2f \n',100*(param(6)-1))
    fprintf('& %6.2f \n',100*Nc/T)
    fprintf('& %6.2f \n',mean(durr(dwo==0)))
    fprintf('& %6.2f \n',mean(dyptr(dwo==0)))
    fprintf('& %6.2f \n',mean(dy2r(dwo==0)))
    fprintf('& %6.2f \n',mean(dkptr(dwo==0)))
    fprintf('& %6.2f \n',mean(dk2r(dwo==0)))
    fprintf('& %6.2f \n',mean(dhpk2r(dwo==0)))
    fprintf('& %6.2f \n',mean(dkm2r(dwo==0)))
    fprintf('& %6.2f \n',mean(dhpkm2r(dwo==0)))
    fprintf('& %6.2f \n',mean(devk2r(dwo==0)))
end
%%

hpy_wc=[
    -.3761499
    -.2445465
    -.3047085
    -.1110554
    1.062039
    2.543716
    -.1820715
    -1.592462
    -.9109716
    -.6551871
    -.3320088
    -.115642
    .1891178
    ];

hpy_woc=[
    -.5219812
    -.2453172
    -.2652525
    -.3011536
    .6215964
    1.999578
    -1.468515
    -.810181
    -.0085641
    -.0776712
    .2178919
    -.1100548
    -.1027005
    ];

hpl_wc=[
    -1.04852
    -1.103903
    -2.18129
    -0.8619354
    3.54215
    3.248002
    1.056448
    -0.3207675
    -0.6601641
    -1.14167
    -0.0362367
    -0.2210597
    -0.7373605
    ];

hpl_woc=[
    -0.7947602
    -1.694641
    -0.4006178
    0.7941448
    1.179271
    0.6100997
    -0.0965231
    -0.3295644
    -0.1883499
    0.1565885
    -0.1126513
    -0.2508744
    -0.4586314
    ];


ylimy=[-6.5 6.5];
yliml=[-6.5 6.5];



if gr_
    ms=4;
    fs=12;
    
    
    subplot(221);h=recplotband((-6:6)',[mean(hpyr(dwo==0,:))' hpy_wc zeros(13,1)],mean(hpyr(dwo==0,:))'-std(hpyr(dwo==0,:))',mean(hpyr(dwo==0,:))'+std(hpyr(dwo==0,:))');
    set(h(1),'linewidth',2,'color','k','marker','o','markeredgecolor','k','markerfacecolor','k','markersize',ms)
    set(h(2),'linewidth',2,'color','k','marker','d','markeredgecolor','k','markerfacecolor','k','markersize',ms,'linestyle','--')
    set(h(3),'linewidth',1,'color','k','linestyle','--')
    set(gca,'fontname','times','fontsize',fs,'xlim',[-6 6],'ylim',ylimy);
    title({'Output','(% deviation about trend)'},'fontname','times','fontsize',fs);
    grid
    
    subplot(222);h=recplotband((-6:6)',[mean(hpkr(dwo==0,:))' hpl_wc zeros(13,1)],mean(hpkr(dwo==0,:))'-std(hpkr(dwo==0,:))',mean(hpkr(dwo==0,:))'+std(hpkr(dwo==0,:))');
    set(h(1),'linewidth',2,'color','k','marker','o','markeredgecolor','k','markerfacecolor','k','markersize',ms)
    set(h(2),'linewidth',2,'color','k','marker','d','markeredgecolor','k','markerfacecolor','k','markersize',ms,'linestyle','--')
    set(h(3),'linewidth',1,'color','k','linestyle','--')
    set(gca,'fontname','times','fontsize',fs,'xlim',[-6 6],'ylim',yliml);
    title({'Credit','(% deviation about trend)'},'fontname','times','fontsize',fs);
    grid
    %     pause
    print('-depsc2','figures/model_vs_data_financial');
    close
    
    subplot(221);h=recplotband((-6:6)',[mean(hpyr(dwo==1,:))' hpy_woc zeros(13,1)],mean(hpyr(dwo==1,:))'-std(hpyr(dwo==1,:))',mean(hpyr(dwo==1,:))'+std(hpyr(dwo==1,:))');
    set(h(1),'linewidth',2,'color','k','marker','o','markeredgecolor','k','markerfacecolor','k','markersize',ms)
    set(h(2),'linewidth',2,'color','k','marker','d','markeredgecolor','k','markerfacecolor','k','markersize',ms,'linestyle','--')
    set(h(3),'linewidth',1,'color','k','linestyle','--')
    set(gca,'fontname','times','fontsize',fs,'xlim',[-6 6],'ylim',ylimy);
    title({'Output','(% deviation about trend)'},'fontname','times','fontsize',fs);
    grid
    
    subplot(222);h=recplotband((-6:6)',[mean(hpkr(dwo==1,:))' hpl_woc zeros(13,1)],mean(hpkr(dwo==1,:))'-std(hpkr(dwo==1,:))',mean(hpkr(dwo==1,:))'+std(hpkr(dwo==1,:))');
    set(h(1),'linewidth',2,'color','k','marker','o','markeredgecolor','k','markerfacecolor','k','markersize',ms)
    set(h(2),'linewidth',2,'color','k','marker','d','markeredgecolor','k','markerfacecolor','k','markersize',ms,'linestyle','--')
    set(h(3),'linewidth',1,'color','k','linestyle','--')
    set(gca,'fontname','times','fontsize',fs,'xlim',[-6 6],'ylim',yliml);
    title({'Credit','(% deviation about trend)'},'fontname','times','fontsize',fs);
    grid
    %     pause
    print('-depsc2','figures/model_vs_data_normal');
    close
    
end