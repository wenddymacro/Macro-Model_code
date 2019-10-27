
function z=soe1_ff(y)
z=zeros(19,1);
global ex_ ex_det_ it_ recur_

global betaus1 betaus2 betaus3 gammaus1 gammaus2 gammaus4 growthus_ss  ...
lambdaus1 lambdaus2 rhous 
global rrus_eq_ss tauus wappaus1 wappaus2 
z(1) = y(32) -(wappaus1*y(13)+wappaus2*y(33)+ex_(it_-3,8));
z(2) = y(32) -(y(31)-y(30));
z(3) = y(31) -(y(12)+ex_(it_-3,7));
z(4) = y(21) -(y(22)+y(33));
z(5) = y(19) -(tauus*growthus_ss+(1-tauus)*y(3)+ex_(it_-3,1));
z(6) = y(22) -(y(5)+y(19)/4+ex_(it_-3,2));
z(7) = y(33) -(betaus1*y(14)+betaus2*y(35)-betaus3*(y(9)-y(10))+ ...
ex_(it_-3,9));
z(8) = y(25) -(lambdaus1*y(37)+(1-lambdaus1)*y(8)+lambdaus2*y(14)+ ...
ex_(it_-3,4));
z(9) = y(29) -(gammaus1*y(11)+(1-gammaus1)*(y(28)+y(37)+gammaus2*(y(37)- ...
y(36))+gammaus4*y(33))+ex_(it_-3,6));
z(10) = y(27) -(y(29)-y(34));
z(11) = y(20) -(y(4)+y(25)/4);
z(12) = y(28) -(rhous*rrus_eq_ss+(1-rhous)*y(10)+ex_(it_-3,5));
z(13) = y(26) -((y(25)+y(7)+y(2)+y(1))/4);
z(14) = y(23) -((y(38)+y(39)+y(40)+y(41)+y(42))/5);
z(15) = y(24) -(y(6)+ex_(it_-3,3));
z(16) = y(18) -(y(37));
z(17) = y(15) -(y(34));
z(18) = y(16) -(y(35));
z(19) = y(17) -(y(36));
