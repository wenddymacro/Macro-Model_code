% [fy, Fall, Udecomp, ndiffuse] = calc_fcast(nahead)
%
% This calculates n-step ahead forecasts for n=1..nahead, storing the
% forecasts (of observational variables) in fy and theoretical n-ahead
% forecast errors in Fall. It returns also Udecomp being a shock
% decomposition of n-ahead forecast update from forecast based on
% information at t-1 to forecast based on information at t.
%
% The 3-dimensional arrays fy and Fall are organized as follows:
%
% fy(:,h,1) is h-ahead forecast based only on initial conditions
% fy(:,h,t) is h-ahead forecast based on information at time t-1
% fy(:,h,end) is h-ahead out-of-sample forecast based on all information
%
% Thus, if gend is a number of observations, then size(fy,3)=gend+1
%
% h-ahead forecast error is then equal to Y(:,t)-fy(:,h,t-h+1)
% 
% The Udecomp is 4-dimensional array and is organized as follows:
% 
% Udecomp(i,h,t,:) is a vector of non-negative fractions summing up to 1
% corresponding to decomposition of change in forecast of Y(i,t+h-1)
% based on information on t-1 and t. The vector corresponds to
% contributions of the shocks at time t-1 (smoothed at time t).
% 
% size(Udecomp,3)=gend
%
% The returned value of ndiffuse says for how many periods there were
% diffuse state estimates. For these periods no forecast is calculated. The
% first meaningful forecast is based on information of ndiffuse
% periods. This means that the first usable value in fy is
% fy(:,:,ndiffuse+1)
%

function [fy, Fall, Udecomp, ndiffuse] = calc_fcast(nahead)

global bayestopt_ exo_nbr dr_ estim_params_ Sigma_e_ options_
global trend_coeff_

nobs = size(options_.varobs,1);
gend = options_.nobs(1);

% preallocate fy and Fall
fy = zeros(nobs, nahead, gend+1);
fx = zeros(size(dr_.ghx,1), nahead, gend+1);
Fall = zeros(nobs, nahead, gend+1);
% preallocate Udecomp
Udecomp = zeros(nobs, nahead, gend, exo_nbr);

% load data
rawdata = read_variables(options_.datafile,options_.varobs,[]);
% transform the data
if options_.loglinear == 1
  rawdata = log(rawdata);
end
if ~isreal(rawdata)
  error(['There are complex values in the data. Probably  a wrong' ...
     ' transformation'])
end
% pickup required data
if options_.prefilter == 1
  bayestopt_.mean_varobs = mean(rawdata(options_.first_obs+(0:gend-1),:),1);
  data = rawdata(options_.first_obs+(0:gend-1),:)'-...
         repmat(bayestopt_.mean_varobs',1,gend);
else
  data = rawdata(options_.first_obs+(0:gend-1),:)';
end  


% set Q, covariance of errors in transition equation
Q = Sigma_e_;
% structiral parameters are supposed to be in the namespace

% solve the model and form state space
[T,R,SteadyState,info] = dynare_resolve;
if info(1) ~= 0
  error('Could not solve the model');
end

% set trend
if options_.loglinear == 1
  constant = log(SteadyState(bayestopt_.mfys));
else
  constant = SteadyState(bayestopt_.mfys);
end
if bayestopt_.with_trend == 1
  trend_coeff = zeros(nobs,1);
  for i=1:nobs
    trend_coeff(i) = evalin('base',bayestopt_.trend_coeff{i});
  end
  trend = constant*ones(1,gend+nahead)+trend_coeff*(1:gend+nahead);
else
  trend = constant*ones(1,gend+nahead);
end

start = options_.presample+1;
np    = size(T,1);
mf    = bayestopt_.mf;

% initial conditions for the kalman filer
if options_.lik_init == 1		% Kalman filter
  Pstar = lyapunov_symm(T,R*Q*transpose(R));
  Pinf	= [];
elseif options_.lik_init == 2	% Old Diffuse Kalman filter
  Pstar = 10*eye(np);
  Pinf	= [];
elseif options_.lik_init == 3	% Diffuse Kalman filter
  Pstar = zeros(np,np);
  ivs = bayestopt_.i_T_var_stable;
  Pstar(ivs,ivs) = lyapunov_symm(T(ivs,ivs),R(ivs,:)*Q* ...
                                 transpose(R(ivs,:)));
  Pinf  = bayestopt_.Pinf;
end

Y = data;

% diffuse kalman filter with forecasts (code of the filter taken almost
% verbatim from DiffuseLikelihood3.m)
mf = bayestopt_.mf;
pp     = size(Y,1);
mm     = size(T,1);
smpl   = size(Y,2);
a      = zeros(mm,1);
QQ     = R*Q*transpose(R);
Ztmp = zeros(mm,pp); Ztmp(mf+mm*[0:pp-1])=1;
Z      = Ztmp';
ZRQinv = inv(Z*QQ*Z');
t      = 0;
crit      	= options_.kalman_tol;
crit1      	= 1.e-6;
newRank	 	= rank(Pinf,crit1);
icc=0;
while newRank & t < smpl
  t = t+1;
  for i=1:pp
    v(i) 	= Y(i,t)-a(mf(i))-trend(i,t);
    Fstar 	= Pstar(mf(i),mf(i));
    Finf	= Pinf(mf(i),mf(i));
    Kstar 	= Pstar(:,mf(i));
    if Finf > crit & newRank,  %added newRank criterion 
      icc=icc+1;
      Kinf	= Pinf(:,mf(i));
      a		= a + Kinf*v(i)/Finf;
      Pstar	= Pstar + Kinf*transpose(Kinf)*Fstar/(Finf*Finf) - ...
	  (Kstar*transpose(Kinf)+Kinf*transpose(Kstar))/Finf;
      Pinf	= Pinf - Kinf*transpose(Kinf)/Finf;
      F(i,t)    = Finf;
      lFF(i,t)  = 0;
      % start new termination criterion for DKF
      if ~isempty(options_.diffuse_d),  
	newRank = (icc<options_.diffuse_d);  
	if newRank & (any(diag(Pinf(mf,mf))>crit)==0 & rank(Pinf,crit1)==0); 
	  options_.diffuse_d = icc;
	  newRank=0;
	end
      else
	newRank = (any(diag(Pinf(mf,mf))>crit) | rank(Pinf,crit1));                 
	if newRank==0, 
	  P0=	T*Pinf*transpose(T);
	  newRank = (any(diag(Pinf(mf,mf))>crit) | rank(P0,crit1));   
	  if newRank==0, 
	    options_.diffuse_d = icc;
	  end
	end                    
      end,
      % end new termination and checks for DKF and fmax
    elseif Fstar > crit 
      %% Note that : (1) rank(Pinf)=0 implies that Finf = 0, (2) outside this loop (when for some i and t the condition
      %% rank(Pinf)=0 is satisfied we have P = Pstar and F = Fstar and (3) Finf = 0 does not imply that
      %% rank(Pinf)=0. [stéphane,11-03-2004].	  
      %if rank(Pinf,crit) == 0
      % the likelihood terms should alwasy be cumulated, not only
      % when Pinf=0, otherwise the lik would depend on the ordering
      % of observed variables
      %end
      F(i,t) = Fstar;
      llF(i,t) = v(i)*v(i)/Fstar;
      a	= a + Kstar*v(i)/Fstar;
      Pstar = Pstar - Kstar*transpose(Kstar)/Fstar;
    else
      % disp(['zero F term in DKF for observed ',int2str(i),' ',num2str(Fi)])
    end
  end 

  if newRank,
    oldRank = rank(Pinf,crit1);
  else
    oldRank = 0;
  end
  a 		= T*a;
  Pstar 	= T*Pstar*transpose(T)+QQ;
  Pinf	= T*Pinf*transpose(T);
  if newRank,
    newRank = rank(Pinf,crit1);  % new crit1 is used 
  end
  if oldRank ~= newRank
    disp('DiffuseLiklihood3 :: T does influence the rank of Pinf!')	
  end  
end
if t == smpl                                                           
  error(['There isn''t enough information to estimate the initial' ... 
	 ' conditions of the nonstationary variables']);                   
end   

% set ndiffuse
ndiffuse = t;

while t < smpl
  t = t+1;
  % now vector a is E(alpha_t|t-1)
  % one step ahead forecast from t-1 will be a(mf) + trend(:,t)
  % h step ahead forecast from t-1 will be (T^{h-1}*a)(mf) + trend(:,t+h-1)
  af = a;
  Pstar1 = Pstar;
  for h=1:nahead
    fy(:,h,t) = af(mf) + trend(:,t+h-1);
    fx(:,h,t) = af(1:size(fx,1));
    af = T*af;
    Fall(:,h,t) = diag(Pstar1(mf,mf));
    Pstar1 = T*Pstar1*transpose(T) + QQ;
  end

  % filter
  v = Y(:,t) - Z*a - trend(:,t);
  Fstar = Z*Pstar*Z';
  K = T*Pstar*Z'/Fstar;
  L = T - K*Z;
  % smooth one step back to get r_tm1t
  r_tm1t = zeros(mm,1);
  r_tm1t = Z'*(Fstar\v) + L'*r_tm1t;
  
  % calculate eta_tm1t
  eta_tm1t = Q*R'*r_tm1t;
  % calculate decomposition
  Ttok = eye(mm,mm); 
  for h = 1:nahead
    for j=1:exo_nbr
      eta=zeros(exo_nbr,1);
      eta(j) = eta_tm1t(j);
      Udecomp(:,h,t,j) = Z*Ttok*Pstar*Z'*ZRQinv*Z*R*eta;
    end
    Ttok = Ttok*T;
  end
  
  a = T*a + K*v;
  Pstar = T*Pstar*L' + QQ;
end

% and calculate forecasts conditional all observations, these are stored
% to fy(:,h,smpl+1)
t = t+1;
af = a;
Pstar1 = Pstar;
for h=1:nahead
  fy(:,h,t) = af(mf) + trend(:,t+h-1);
  fx(:,h,t) = af(1:size(fx,1));
  af = T*af;
  Fall(:,h,t) = diag(Pstar1(mf,mf));
  Pstar1 = T*Pstar1*transpose(T) + QQ;
end
