% ** ONLY UNDER UNIX SYSTEM
%path(path,'/usr2/f1taz14/mymatlab')

if (1)
   close all
   clear all
   clc  
end

%-- path for restriction .m files
%path(path,'D:\ZhaData\TZCcode\Model_TVSVAR\Identification'); %For Windows.
%path(path,'/Users/tzha/ZhaData/TZCcode/Model_MSBVAR/Identification');


%===========================================
% Exordium I
%===========================================
format short g     % format
%

%*** Load data and series

% data_in = xlsread('CLSZ_data.xlsx', 'VAR_data', 'B2:Z150');  
% 
% save('clsz_data_in', 'data_in');
% clear ('data_in');
% return;
               
% data_in = xlsread('zha_data.xlsx', 'for_BVAR', 'B2:Z200');  % 1990:Q1 - 2015:Q4
% 
% save('zha_data_in', 'data_in');
% clear ('data_in');
% % 
% return;
% % 


%------ Variables in the data loaded (zha_data_in.mat) --------
%
% 1: RRR	
% 2: Benchmark lending rate (one year)
% 3: Heavy industry loans in log units (2003Q1-2015Q4)
% 4: Light industry loans in log units (2003Q1 - 2015Q4)
% 5: Share of heavey industry loan in total loans (2003Q1- 2015Q4)
% 6: Nominal SOE GCF (1995Q1-2015Q4)
% 7. Nominal Business GCF = Nominal SOE GCF +  Nominal
%   Private GCF + Nominal non-SOE GCF  (1995Q1-2015Q4)
% 8. SOE investment share = (6)/(7), 1995Q1-2015Q4


% ----------- Prompt selection of a shock and set up sample ranges --------

indx_model = 2;
load clsz_data_in;
if indx_model == 1;
    q_m = 4;   % quarters (4) or months (12)
    yrBin=2003;   % beginning of the year
    qmBin=1;    % begining of the quarter or month
    yrFin=2015;   % final year
    qmFin=4;    % final month or quarter
    nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
    xdd = data_in(73:124,[1 3 8 2]);  
    indx_numvars = 4;  
    vlist=[1:4];
    varlist={'Required reserve ratio', 'Heavy industry loan (log)', 'SOE investment share', 'Lending rate'};
    vlistlog = [];   % subset of "vlist.  Variables in log level so that differences are in *quraterly* growth, unlike L which is in level.
    vlistper = [1:4];     % subset of "vlist"
    rnum = 2;      % number of rows in the graph
    cnum = 2;      % number of columns in the graph
elseif indx_model == 2;
    q_m = 4;   % quarters (4) or months (12)
    yrBin=2003;   % beginning of the year
    qmBin=1;    % begining of the quarter or month
    yrFin=2015;   % final year
    qmFin=4;    % final month or quarter
    nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
    xdd = data_in(73:124,[3 8 2 1]);  
    indx_numvars = 4;  
    vlist=[1:4];
    varlist={'Heavy industry loan (log)', 'SOE investment share', 'Lending rate', 'Required reserve ratio'};
    vlistlog = [1];   % subset of "vlist.  Variables in log level so that differences are in *quraterly* growth, unlike L which is in level.
    vlistper = [2:4];     % subset of "vlist"
    rnum = 2;      % number of rows in the graph
    cnum = 2;      % number of columns in the graph
elseif indx_model == 3;
    q_m = 4;   % quarters (4) or months (12)
    yrBin=2003;   % beginning of the year
    qmBin=1;    % begining of the quarter or month
    yrFin=2015;   % final year
    qmFin=4;    % final month or quarter
    nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
    xdd = data_in(73:124,[1 2 5]);  
    indx_numvars = 3;  
    vlist=[1:3];
    varlist={'Required reserve ratio', 'Lending rate', 'Heavy industry loan share'};
    vlistlog = [];   % subset of "vlist.  Variables in log level so that differences are in *quraterly* growth, unlike L which is in level.
    vlistper = [1:3];     % subset of "vlist"
    rnum = 2;      % number of rows in the graph
    cnum = 2;      % number of columns in the graph
elseif indx_model == 4;
    q_m = 4;   % quarters (4) or months (12)
    yrBin=2003;   % beginning of the year
    qmBin=1;    % begining of the quarter or month
    yrFin=2015;   % final year
    qmFin=4;    % final month or quarter
    nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
    xdd = data_in(73:124,[5 2 1]);  
    indx_numvars = 3;  
    vlist=[1:3];
    varlist={'Heavy industry loan share', 'Lending rate', 'Required reserve ratio'};
    vlistlog = [];   % subset of "vlist.  Variables in log level so that differences are in *quraterly* growth, unlike L which is in level.
    vlistper = [1:3];     % subset of "vlist"
    rnum = 2;      % number of rows in the graph
    cnum = 2;      % number of columns in the graph
elseif indx_model==5;
    q_m = 4;   % quarters (4) or months (12)
    yrBin=2003;   % beginning of the year
    qmBin=1;    % begining of the quarter or month
    yrFin=2015;   % final year
    qmFin=4;    % final month or quarter
    nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
    xdd = data_in(73:124,[1 2 3 8]);  
    indx_numvars = 4;  
    vlist=[1:4];
    varlist={'Required reserve ratio', 'Lending rate', 'Heavy industry loan (log)', 'SOE investment share'};
    vlistlog = [3];   % subset of "vlist.  Variables in log level so that differences are in *quraterly* growth, unlike L which is in level.
    vlistper = [1 2 4];     % subset of "vlist"
    rnum = 2;      % number of rows in the graph
    cnum = 2;      % number of columns in the graph
elseif indx_model==6;
    q_m = 4;   % quarters (4) or months (12)
    yrBin=1999;   % beginning of the year
    qmBin=1;    % begining of the quarter or month
    yrFin=2015;   % final year
    qmFin=4;    % final month or quarter
    nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
    xdd = data_in(37:104,[1,10,14]);  
    indx_numvars = 3;  
    vlist=[1:3];
    varlist={'Required reserve ratio', 'SOE revenue share', 'Interest rate'};
    vlistlog = [];   % subset of "vlist.  Variables in log level so that differences are in *quraterly* growth, unlike L which is in level.
    vlistper = [1:3];     % subset of "vlist"
    rnum = 2;      % number of rows in the graph
    cnum = 2;      % number of columns in the graph
elseif indx_model == 9;
    q_m = 4;   % quarters (4) or months (12)
    yrBin=1995;   % beginning of the year
    qmBin=1;    % begining of the quarter or month
    yrFin=2013;   % final year
    qmFin=4;    % final month or quarter
    nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
    xdd = data_in(21:96,[1,14,5,16]);  
    indx_numvars = 4;  
    vlist=[1:4];
    varlist={'Required reserve ratio', 'Interest rate', 'GDP','SOE investment share'};
    vlistlog = [];   % subset of "vlist.  Variables in log level so that differences are in *quraterly* growth, unlike L which is in level.
    vlistper = [1:4];     % subset of "vlist"
    rnum = 2;      % number of rows in the graph
    cnum = 2;      % number of columns in the graph
elseif indx_model == 10;
    q_m = 4;   % quarters (4) or months (12)
    yrBin=1995;   % beginning of the year
    qmBin=1;    % begining of the quarter or month
    yrFin=2013;   % final year
    qmFin=4;    % final month or quarter
    nData=(yrFin-yrBin)*q_m + (qmFin-qmBin+1);
    xdd = data_in(21:96,[14,16,5,1]);  
    indx_numvars = 4;  
    vlist=[1:4];
    varlist={'Interest rate', 'GDP','SOE investment share', 'Required reserve ratio'};
    vlistlog = [];   % subset of "vlist.  Variables in log level so that differences are in *quraterly* growth, unlike L which is in level.
    vlistper = [1:4];     % subset of "vlist"
    rnum = 2;      % number of rows in the graph
    cnum = 2;      % number of columns in the graph
end;


[nt,ndv]=size(xdd);
if nt~=nData
   disp(' ')
   warning(sprintf('nt=%d, Caution: not equal to the length in the data',nt));
   %disp(sprintf('nt=%d, Caution: not equal to the length in the data',nt));
   disp('Press ctrl-c to abort')
   return
end   

idfile_const='fn_iden_upperchol';   %Only used by msstart2.m.

% idfile_const='fn_iden_lowerchol';   %Only used by msstart2.m.
%
% rnum = 2;      % number of rows in the graph
% cnum = 2;      % number of columns in the graph
ylab = varlist;
xlab = varlist;



% indx_1lead2 = 1; %1: Housing price (1) ordered ahead of unemployment (2); 0: the other way around.
% if (indx_1lead2)
%   vlist = [1 2];
%   varlist={'Housing price', 'Unemployment'};
%   xlab = {'Shock to housing price','Shock to unemployment'};      
% else   
%   vlist = [2 1];
%   varlist={'Unemployment','Housing price'};                       
%   xlab = {'Shock to unemployment','Shock to housing price'};      
% end   
% vlistlog = [1:2];   % subset of "vlist.  Variables in log level so that differences are in *quraterly* growth, unlike L which is in level.
% vlistper = [];           % subset of "vlist"
% idfile_const='ftd_upperchol2v';   %Only used by msstart2.m.
%    %
% ylab = varlist;
% rnum = 2;      % number of rows in the graph
% cnum = 1;      % number of columns in the graph

%----------------
nvar = length(vlist);   % number of endogenous variables
nlogeno = length(vlistlog)  % number of endogenous variables in vlistlog
npereno = length(vlistper)  % number of endogenous variables in vlistper
if (nvar~=(nlogeno+npereno))
   disp(' ')
   warning('Check xlab, nlogeno or npereno to make sure of endogenous variables in vlist')
   disp('Press ctrl-c to abort')
   return
elseif (nvar==length(vlist))
   nexo=1;    % only constants as an exogenous variable.  The default setting.
elseif (nvar<length(vlist))
   nexo=length(vlist)-nvar+1;
else
   disp(' ')
   warning('Make sure there are only nvar endogenous variables in vlist')
   disp('Press ctrl-c to abort')
   return
end
if rnum*cnum<nvar
   warning('rnum*cnum must be at least as large as nvar')
   disp('Hit any key to continue, or ctrl-c to abort')
   pause
end


%------- A specific sample is considered for estimation -------
yrStart=yrBin;
qmStart=qmBin;
yrEnd=yrFin; 
qmEnd=qmFin; %1;
qm_dates_all_sample = ((yrStart+(qmStart-1)/q_m):1/q_m:(yrEnd+(qmEnd-1)/q_m))';
nfyr = 4;   % number of years for forecasting
if nfyr<1
   error('To be safe, the number of forecast years should be at least 1')
end
ystr=num2str(yrEnd);
forelabel = [ ystr(3:4) ':' num2str(qmEnd) ' Forecast'];

nSample=(yrEnd-yrStart)*q_m + (qmEnd-qmStart+1);
if qmEnd==q_m
   E1yrqm = [yrEnd+1 1];  % first year and quarter (month) after the sample
else
   E1yrqm = [yrEnd qmEnd+1];  % first year and quarter (month) after the sample
end
E2yrqm = [yrEnd+nfyr qmEnd];   % end at the last month (quarter) of a calendar year after the sample
[fdates,nfqm]=fn_calyrqm(q_m,E1yrqm,E2yrqm);   % forecast dates and number of forecast dates
[sdates,nsqm] = fn_calyrqm(q_m,[yrStart qmStart],[yrEnd qmEnd]);
   % sdates: dates for the whole sample (including lags)
if nSample~=nsqm
   warning('Make sure that nSample is consistent with the size of sdates')
   disp('Hit any key to continue, or ctrl-c to abort')
   pause
end
imstp = 4*q_m;    % <<>>  impulse responses (4 years)
nayr = 4; %nfyr;  % number of years before forecasting for plotting.


%------- Prior, etc. -------
lags = q_m;        % number of lags
qm_dates_ess = qm_dates_all_sample(lags+1:end);
indxC0Pres = 0;   % 1: cross-A0-and-A+ restrictions; 0: idfile_const is all we have
            % Example for indxOres==1: restrictions of the form P(t) = P(t-1).
Rform = 0;  % 1: contemporaneous recursive reduced form; 0: restricted (non-recursive) form
Pseudo = 1;  % 1: Pseudo forecasts; 0: real time forecasts
indxPrior = 0;  % 1: Bayesian prior; 0: no prior
indxDummy = indxPrior;  % 1: add dummy observations to the data; 0: no dummy added.
ndobs = 0;  % No dummy observations for xtx, phi, fss, xdatae, etc.  Dummy observations are used as an explicit prior in fn_rnrprior_covres_dobs.m.
%if indxDummy
%   ndobs=nvar+1;         % number of dummy observations
%else
%   ndobs=0;    % no dummy observations
%end
%=== The following mu is effective only if indxPrior==1.
mu = zeros(6,1);   % hyperparameters
mu(1) = 1;
mu(2) = 1;
mu(3) = 1;
mu(4) = 1.2;
mu(5) = 1;
mu(6) = 1;
%mu(1) = 1;
%mu(2) = 0.5;
%mu(3) = 0.1;
%mu(4) = 1.0;
%mu(5) = 1.0;
%mu(6) = 1.0;
%   mu(1): overall tightness and also for A0;
%   mu(2): relative tightness for A+;
%   mu(3): relative tightness for the constant term;
%   mu(4): tightness on lag decay;  (1)
%   mu(5): weight on nvar sums of coeffs dummy observations (unit roots);
%   mu(6): weight on single dummy initial observation including constant
%           (cointegration, unit roots, and stationarity);
%
%
hpmsmd = [0.0; 0.0];
indxmsmdeqn = [1; 2; 1; 2];


tdf = 3;          % degrees of freedom for t-dist for initial draw of the MC loop
nbuffer = 200;        % a block or buffer of draws (buffer) that is saved to the disk (not memory)
ndraws1=1*nbuffer;         % 1st part of Monte Carlo draws
ndraws2=10*ndraws1         % 2nd part of Monte Carlo draws
seednumber = 0; %7910;    %472534;   % if 0, random state at each clock time
           % good one 420 for [29 45], [29 54]
if seednumber
   randn('state',seednumber);
   rand('state',seednumber);
else
   randn('state',fix(100*sum(clock)));
   rand('state',fix(100*sum(clock)));
end
%  nstarts=1         % number of starting points
%  imndraws = nstarts*ndraws2;   % total draws for impulse responses or forecasts
%<<<<<<<<<<<<<<<<<<<




