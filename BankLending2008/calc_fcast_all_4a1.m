% adapted for Dynare version 4. February 2008, Michel Juillard
%
% [fy, fx, Fyall, Fxall, Udecomp, ndiffuse] = calc_fcast_al_4()
%
% This calculates n-step ahead forecasts for n=1..nahead, storing the
% forecasts (of observational variables) in fy and theoretical n-ahead
% forecast errors in Fall. It returns also Udecomp being a shock
% decomposition of n-ahead forecast update from forecast based on
% information at t-1 to forecast based on information at t.
%
% **** fy ****
% The 3-dimensional array fy is organized as follows:
%
% the first dimension of fy corresponds to observed variables in the
%           options_.varobs
% fy(:,h,1) is h-ahead forecast based only on initial conditions
% fy(:,h,t) is h-ahead forecast based on information at time t-1
% fy(:,h,end) is h-ahead out-of-sample forecast based on all information
%
% Thus, if gend is a number of observations, then size(fy,3)=gend+1
%
% h-ahead forecast error is then equal to Y(:,t)-fy(:,h,t-h+1)
%
% **** fx ****
% The 3-dimensional array fx is organized in the same way as fx, the only
% difference is that it contains all variables in the dr_.order_var. So
% the ordering of the variables (first dimension) can be displayed as
% lgy_(dr_.order_var). The values do not have the steady state and
% trends added.
%
% **** Fyall ****
% The 3-dimensional array Fyall is organized as fy, it contains theoretical
% stderrors of the forecasts
%
% **** Fxall ****
% The 3-dimensianal array Fxall is organized as fx, it contains
% theoretical stderrors pf the forecasts
%
% **** Udecomp ****
% The Udecomp is 4-dimensional array and is organized as follows:
%
% Udecomp(i,h,t,:) is a vector of non-negative fractions summing up to 1
% corresponding to decomposition of change in forecast of Y(i,t+h-1)
% based on information on t-1 and t. The vector corresponds to
% contributions of the shocks at time t-1 (smoothed at time t).
%
% size(Udecomp,3)=gend
%
% **** ndiffuse ****
% The returned value of ndiffuse says for how many periods there were
% diffuse state estimates. For these periods no forecast is calculated. The
% first meaningful forecast is based on information of ndiffuse
% periods. This means that the first usable value in fy is
% fy(:,:,ndiffuse+1)
%

function [fy, fx, Fyall, Fxall, Udecomp, ndiffuse] = calc_fcast_all_1()

global options_ oo_ bayestopt_ M_


nobs = size(options_.varobs,1);
gend = options_.nobs(1);
nk = length(options_.filter_step_ahead);
exo_nbr = M_.exo_nbr;

% preallocate fy and Fall
fy = zeros(nobs, nk, gend+1);
fx = zeros(size(oo_.dr.ghx,1), nk, gend+1);
Fyall = zeros(nobs, nk, gend+1);
Fxall = zeros(size(oo_.dr.ghx,1), nk, gend+1);

% preallocate Udecomp
Udecomp = zeros(size(oo_.dr.ghx,1), nk, gend, exo_nbr);

mf = bayestopt_.mf;
mfys = bayestopt_.mfys;

trend = repmat(oo_.Smoother.SteadyState(mfys),1,gend+nk);

if bayestopt_.with_trend
   trend = trend + oo_.Smoother.TrendCoeffs*[1:gend+nk];
end

for t=2:gend+1
   for jnk = 1:nk
       fx(:,jnk,t) = oo_.FilteredVariablesKStepAhead(jnk,:,t+jnk-1);
       fy(:,jnk,t) = squeeze(oo_.FilteredVariablesKStepAhead(jnk,mf,t+jnk-1))'+trend(:,t+jnk-1);
       Fxall(:,jnk,t) = diag(squeeze(oo_.FilteredVariablesKStepAheadVariances(jnk,:,:,t+jnk-1)));
       Fyall(:,jnk,t) = diag(squeeze(oo_.FilteredVariablesKStepAheadVariances(jnk,mf,mf,t+ ...
                                                 jnk-1)));
       if t <= gend
           Udecomp(:,jnk,t,:) = ...
               oo_.FilteredVariablesShockDecomposition(jnk,:,:,t+jnk-1);
       end
   end
end

ndiffuse = oo_.Smoother.integration_order;