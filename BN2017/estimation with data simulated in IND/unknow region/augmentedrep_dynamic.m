function [residual, g1, g2, g3] = augmentedrep_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(10, 1);
beta__ = (1+params(8)/100)^(-0.25);
tau__ = 1/params(1);
lhs =y(6);
rhs =y(16)-tau__*(y(7)-y(17))+y(9);
residual(1)= lhs-rhs;
lhs =y(8);
rhs =y(17)*beta__+params(2)*(y(6)-y(10));
residual(2)= lhs-rhs;
lhs =y(7);
rhs =params(3)*y(1)+(1-params(3))*(y(8)*params(6)+(y(6)-y(10))*params(7))+x(it_, 1);
residual(3)= lhs-rhs;
lhs =y(9);
rhs =params(4)*y(2)+x(it_, 2);
residual(4)= lhs-rhs;
lhs =y(10);
rhs =params(5)*y(3)+x(it_, 3);
residual(5)= lhs-rhs;
lhs =y(14);
rhs =y(17);
residual(6)= lhs-rhs;
lhs =y(15);
rhs =params(10)*y(5)+x(it_, 4)-(y(8)-y(4));
residual(7)= lhs-rhs;
lhs =y(11);
rhs =y(6);
residual(8)= lhs-rhs;
lhs =y(13);
rhs =params(9)+4*y(8);
residual(9)= lhs-rhs;
lhs =y(12);
rhs =params(8)+params(9)+4*y(7);
residual(10)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(10, 21);

  %
  % Jacobian matrix
  %

  g1(1,6)=1;
  g1(1,16)=(-1);
  g1(1,7)=tau__;
  g1(1,17)=(-tau__);
  g1(1,9)=(-1);
  g1(2,6)=(-params(2));
  g1(2,8)=1;
  g1(2,17)=(-beta__);
  g1(2,10)=params(2);
  g1(3,6)=(-((1-params(3))*params(7)));
  g1(3,1)=(-params(3));
  g1(3,7)=1;
  g1(3,8)=(-((1-params(3))*params(6)));
  g1(3,10)=(-((1-params(3))*(-params(7))));
  g1(3,18)=(-1);
  g1(4,2)=(-params(4));
  g1(4,9)=1;
  g1(4,19)=(-1);
  g1(5,3)=(-params(5));
  g1(5,10)=1;
  g1(5,20)=(-1);
  g1(6,17)=(-1);
  g1(6,14)=1;
  g1(7,8)=1;
  g1(7,4)=(-1);
  g1(7,5)=(-params(10));
  g1(7,15)=1;
  g1(7,21)=(-1);
  g1(8,6)=(-1);
  g1(8,11)=1;
  g1(9,8)=(-4);
  g1(9,13)=1;
  g1(10,7)=(-4);
  g1(10,12)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,441);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],10,9261);
end
end
end
end
