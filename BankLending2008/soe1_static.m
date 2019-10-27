function [residual, g1, g2] = soe1_static( y, x )
% 
% Status : Computes static model for Dynare
% 
% Warning : this file is generated automatically by Dynare
%   from model file (.mod)

global M_ 
if M_.param_nbr > 0
  params = M_.params;
end
  residual = zeros( 24, 1);

	% 
	% Model equations
	% 

lhs =y(4);
rhs =params(3)*y(4)+params(4)*y(8)+x(2);
residual(1)= lhs-rhs;
lhs =y(4);
rhs =y(5)-y(3);
residual(2)= lhs-rhs;
lhs =y(5);
rhs =y(5)+y(17)+x(3);
residual(3)= lhs-rhs;
lhs =y(17);
rhs =(1-params(16))*y(17)+x(9);
residual(4)= lhs-rhs;
lhs =y(9);
rhs =y(8)+y(10);
residual(5)= lhs-rhs;
lhs =y(12);
rhs =params(5)*params(6)+(1-params(5))*y(12)+x(5);
residual(6)= lhs-rhs;
lhs =y(10);
rhs =y(10)+y(12)/4+x(7);
residual(7)= lhs-rhs;
lhs =y(8);
rhs =params(7)*y(8)+params(8)*y(8)-params(9)*(y(1)-y(2))-params(18)*(0.04*(y(23)+y(23))+0.08*(y(23)+y(23))+0.12*(y(23)+y(23))+0.16*(y(23)+y(23))+0.2*y(23))+x(6);
residual(8)= lhs-rhs;
lhs =y(23);
rhs =x(10);
residual(9)= lhs-rhs;
lhs =y(21);
rhs =x(10)+y(22)-params(17)*y(8);
residual(10)= lhs-rhs;
lhs =y(22);
rhs =y(22)+x(11);
residual(11)= lhs-rhs;
lhs =y(18);
rhs =4*(y(9)-y(9));
residual(12)= lhs-rhs;
lhs =y(19);
rhs =y(9)-y(9);
residual(13)= lhs-rhs;
lhs =y(20);
rhs =y(10)-y(10);
residual(14)= lhs-rhs;
lhs =y(6);
rhs =params(10)*y(7)+(1-params(10))*y(7)+y(8)*params(11)-x(8);
residual(15)= lhs-rhs;
lhs =y(11);
rhs =params(12)*y(11)+(1-params(12))*(y(2)+y(7)+params(13)*(y(7)-params(15))+y(8)*params(14))+x(4);
residual(16)= lhs-rhs;
lhs =y(1);
rhs =y(11)-y(6);
residual(17)= lhs-rhs;
lhs =y(13);
rhs =y(13)+y(6)/4;
residual(18)= lhs-rhs;
lhs =y(2);
rhs =params(1)*params(2)+y(2)*(1-params(1))+x(1);
residual(19)= lhs-rhs;
lhs =y(7);
rhs =(y(6)+y(6)+y(6)+y(6))/4;
residual(20)= lhs-rhs;
lhs =y(14);
rhs =y(7);
residual(21)= lhs-rhs;
lhs =y(16);
rhs =y(6);
residual(22)= lhs-rhs;
lhs =y(15);
rhs =y(8);
residual(23)= lhs-rhs;
lhs =y(24);
rhs =params(18)*(0.04*(y(23)+y(23))+0.08*(y(23)+y(23))+0.12*(y(23)+y(23))+0.16*(y(23)+y(23))+0.2*y(23));
residual(24)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(24, 24);

	% 
	% Jacobian matrix
	% 

  g1(1, 4)=  g1(1, 4)+1;
  g1(1, 4)=  g1(1, 4)+(-params(3));
  g1(1, 8)=  g1(1, 8)+(-params(4));
  g1(2, 4)=  g1(2, 4)+1;
  g1(2, 5)=  g1(2, 5)+(-1);
  g1(2, 3)=  g1(2, 3)+1;
  g1(3, 5)=  g1(3, 5)+1;
  g1(3, 5)=  g1(3, 5)+(-1);
  g1(3, 17)=  g1(3, 17)+(-1);
  g1(4, 17)=  g1(4, 17)+1;
  g1(4, 17)=  g1(4, 17)+(-(1-params(16)));
  g1(5, 8)=  g1(5, 8)+(-1);
  g1(5, 9)=  g1(5, 9)+1;
  g1(5, 10)=  g1(5, 10)+(-1);
  g1(6, 12)=  g1(6, 12)+1;
  g1(6, 12)=  g1(6, 12)+(-(1-params(5)));
  g1(7, 10)=  g1(7, 10)+1;
  g1(7, 12)=  g1(7, 12)+(-0.25);
  g1(7, 10)=  g1(7, 10)+(-1);
  g1(8, 8)=  g1(8, 8)+1;
  g1(8, 8)=  g1(8, 8)+(-params(7));
  g1(8, 8)=  g1(8, 8)+(-params(8));
  g1(8, 1)=  g1(8, 1)+params(9);
  g1(8, 2)=  g1(8, 2)+(-params(9));
  g1(8, 23)=  g1(8, 23)+params(18)*0.04;
  g1(8, 23)=  g1(8, 23)+params(18)*0.04;
  g1(8, 23)=  g1(8, 23)+params(18)*0.08;
  g1(8, 23)=  g1(8, 23)+params(18)*0.08;
  g1(8, 23)=  g1(8, 23)+params(18)*0.12;
  g1(8, 23)=  g1(8, 23)+params(18)*0.12;
  g1(8, 23)=  g1(8, 23)+params(18)*0.16;
  g1(8, 23)=  g1(8, 23)+params(18)*0.16;
  g1(8, 23)=  g1(8, 23)+params(18)*0.2;
  g1(9, 23)=  g1(9, 23)+1;
  g1(10, 21)=  g1(10, 21)+1;
  g1(10, 22)=  g1(10, 22)+(-1);
  g1(10, 8)=  g1(10, 8)+params(17);
  g1(11, 22)=  g1(11, 22)+1;
  g1(11, 22)=  g1(11, 22)+(-1);
  g1(12, 9)=  g1(12, 9)+(-4);
  g1(12, 18)=  g1(12, 18)+1;
  g1(12, 9)=  g1(12, 9)+4;
  g1(13, 9)=  g1(13, 9)+(-1);
  g1(13, 19)=  g1(13, 19)+1;
  g1(13, 9)=  g1(13, 9)+1;
  g1(14, 10)=  g1(14, 10)+(-1);
  g1(14, 20)=  g1(14, 20)+1;
  g1(14, 10)=  g1(14, 10)+1;
  g1(15, 8)=  g1(15, 8)+(-params(11));
  g1(15, 6)=  g1(15, 6)+1;
  g1(15, 7)=  g1(15, 7)+(-params(10));
  g1(15, 7)=  g1(15, 7)+(-(1-params(10)));
  g1(16, 8)=  g1(16, 8)+(-((1-params(12))*params(14)));
  g1(16, 11)=  g1(16, 11)+1;
  g1(16, 11)=  g1(16, 11)+(-params(12));
  g1(16, 2)=  g1(16, 2)+(-(1-params(12)));
  g1(16, 7)=  g1(16, 7)+(-((1-params(12))*(1+params(13))));
  g1(17, 11)=  g1(17, 11)+(-1);
  g1(17, 1)=  g1(17, 1)+1;
  g1(17, 6)=  g1(17, 6)+1;
  g1(18, 6)=  g1(18, 6)+(-0.25);
  g1(18, 13)=  g1(18, 13)+1;
  g1(18, 13)=  g1(18, 13)+(-1);
  g1(19, 2)=  g1(19, 2)+(-(1-params(1)));
  g1(19, 2)=  g1(19, 2)+1;
  g1(20, 6)=  g1(20, 6)+(-0.25);
  g1(20, 7)=  g1(20, 7)+1;
  g1(20, 6)=  g1(20, 6)+(-0.25);
  g1(20, 6)=  g1(20, 6)+(-0.25);
  g1(20, 6)=  g1(20, 6)+(-0.25);
  g1(21, 7)=  g1(21, 7)+(-1);
  g1(21, 14)=  g1(21, 14)+1;
  g1(22, 6)=  g1(22, 6)+(-1);
  g1(22, 16)=  g1(22, 16)+1;
  g1(23, 8)=  g1(23, 8)+(-1);
  g1(23, 15)=  g1(23, 15)+1;
  g1(24, 23)=  g1(24, 23)+(-(params(18)*0.04));
  g1(24, 23)=  g1(24, 23)+(-(params(18)*0.04));
  g1(24, 23)=  g1(24, 23)+(-(params(18)*0.08));
  g1(24, 23)=  g1(24, 23)+(-(params(18)*0.08));
  g1(24, 23)=  g1(24, 23)+(-(params(18)*0.12));
  g1(24, 23)=  g1(24, 23)+(-(params(18)*0.12));
  g1(24, 23)=  g1(24, 23)+(-(params(18)*0.16));
  g1(24, 23)=  g1(24, 23)+(-(params(18)*0.16));
  g1(24, 23)=  g1(24, 23)+(-(params(18)*0.2));
  g1(24, 24)=  g1(24, 24)+1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
