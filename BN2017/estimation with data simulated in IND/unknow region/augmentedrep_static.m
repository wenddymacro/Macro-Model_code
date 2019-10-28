function [residual, g1, g2, g3] = augmentedrep_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 10, 1);

%
% Model equations
%

beta__ = (1+params(8)/100)^(-0.25);
tau__ = 1/params(1);
lhs =y(1);
rhs =y(1)-tau__*(y(2)-y(3))+y(4);
residual(1)= lhs-rhs;
lhs =y(3);
rhs =y(3)*beta__+params(2)*(y(1)-y(5));
residual(2)= lhs-rhs;
lhs =y(2);
rhs =y(2)*params(3)+(1-params(3))*(y(3)*params(6)+(y(1)-y(5))*params(7))+x(1);
residual(3)= lhs-rhs;
lhs =y(4);
rhs =y(4)*params(4)+x(2);
residual(4)= lhs-rhs;
lhs =y(5);
rhs =y(5)*params(5)+x(3);
residual(5)= lhs-rhs;
lhs =y(9);
rhs =y(3);
residual(6)= lhs-rhs;
lhs =y(10);
rhs =y(10)*params(10)+x(4)-(y(3)-y(9));
residual(7)= lhs-rhs;
lhs =y(6);
rhs =y(1);
residual(8)= lhs-rhs;
lhs =y(8);
rhs =params(9)+y(3)*4;
residual(9)= lhs-rhs;
lhs =y(7);
rhs =params(8)+params(9)+y(2)*4;
residual(10)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(10, 10);

  %
  % Jacobian matrix
  %

  g1(1,2)=tau__;
  g1(1,3)=(-tau__);
  g1(1,4)=(-1);
  g1(2,1)=(-params(2));
  g1(2,3)=1-beta__;
  g1(2,5)=params(2);
  g1(3,1)=(-((1-params(3))*params(7)));
  g1(3,2)=1-params(3);
  g1(3,3)=(-((1-params(3))*params(6)));
  g1(3,5)=(-((1-params(3))*(-params(7))));
  g1(4,4)=1-params(4);
  g1(5,5)=1-params(5);
  g1(6,3)=(-1);
  g1(6,9)=1;
  g1(7,3)=1;
  g1(7,9)=(-1);
  g1(7,10)=1-params(10);
  g1(8,1)=(-1);
  g1(8,6)=1;
  g1(9,3)=(-4);
  g1(9,8)=1;
  g1(10,2)=(-4);
  g1(10,7)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],10,100);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],10,1000);
end
end
end
end
