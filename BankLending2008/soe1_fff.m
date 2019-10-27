
function z=soe1_fff(y)
z=zeros(19,1);
global ex_ ex_det_ it_ recur_

global betaus1 betaus2 betaus3 gammaus1 gammaus2 gammaus4 growthus_ss  ...
lambdaus1 lambdaus2 rhous 
global rrus_eq_ss tauus wappaus1 wappaus2 
z(1) = y(18) -(wappaus1*y(18)+wappaus2*y(19)+ex_(it_-3,8));
z(2) = y(18) -(y(17)-y(16));
z(3) = y(17) -(y(17)+ex_(it_-3,7));
z(4) = y(7) -(y(8)+y(19));
z(5) = y(5) -(tauus*growthus_ss+(1-tauus)*y(5)+ex_(it_-3,1));
z(6) = y(8) -(y(8)+y(5)/4+ex_(it_-3,2));
z(7) = y(19) -(betaus1*y(19)+betaus2*y(19)-betaus3*(y(13)-y(14))+ ...
ex_(it_-3,9));
z(8) = y(11) -(lambdaus1*y(12)+(1-lambdaus1)*y(12)+lambdaus2*y(19)+ ...
ex_(it_-3,4));
z(9) = y(15) -(gammaus1*y(15)+(1-gammaus1)*(y(14)+y(12)+gammaus2*(y(12)- ...
y(10))+gammaus4*y(19))+ex_(it_-3,6));
z(10) = y(13) -(y(15)-y(11));
z(11) = y(6) -(y(6)+y(11)/4);
z(12) = y(14) -(rhous*rrus_eq_ss+(1-rhous)*y(14)+ex_(it_-3,5));
z(13) = y(12) -((y(11)+y(11)+y(11)+y(11))/4);
z(14) = y(9) -((y(12)+y(12)+y(12)+y(12)+y(12))/5);
z(15) = y(10) -(y(10)+ex_(it_-3,3));
z(16) = y(4) -(y(12));
z(17) = y(1) -(y(11));
z(18) = y(2) -(y(19));
z(19) = y(3) -(y(10));
if ~isreal(z)
  z = real(z)+imag(z).^2;
end
