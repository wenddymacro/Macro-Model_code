% Bilbiie, Ghironi, and Melitz (2012), 
%¡°Endogenous Entry, Product Variety, and Business Cycles¡±
%bgm_simp.mod
% Non-capital version of bilbiie, ghironi, melitz (2012)


var rhoo d w l ne n v c z  ;

varexo e;

parameters thetta betta fe chii psii phii deltaf r muu;

betta = .99;
deltaf = .025;
thetta = 3.8;
fe = 1;
chii = .924271;
psii = 4;
phii = .979;

r = betta^(-1)-1;

muu = (thetta)/(thetta-1);
  
model;

% Table 2 equations
% Recall that N is predetermined at time t!

rhoo = muu*w/z;

rhoo = n^(1/(thetta-1));

d = (1-1/muu)*(c/n);

v = w*fe/z;

n = (1-deltaf)*(n(-1) + ne(-1));

chii*l^(1/psii) = w/c;

v = betta*(1-deltaf)*(c/c(+1))*(v(+1)+d(+1));

c + ne*v = w*l + n*d;

% Additional Equations

ln(z) = phii*ln(z(-1)) + e;

% Empirically relevant, log


end;

%initval;
 steady_state_model;

z = 1;

n = ((1-deltaf)/(chii*(thetta)*(r+deltaf)))*((chii*(thetta)*(r+deltaf))/(thetta*(r+deltaf)-r))^(1/(1+psii))*z/fe;

c = (r+deltaf)*(thetta-1)*fe*n^(thetta/(thetta-1))/(1-deltaf);

rhoo = n^(1/(thetta-1));

w = rhoo*z/muu;

d = (1-1/muu)*(c/n);

v = w*fe/z;

ne = (deltaf/(1-deltaf))*n;

l = ((w/c)/chii)^(psii);


end;

steady;
check;

shocks;
var e = .0072^2;
end;

stoch_simul (order=1,irf=5) l rhoo  ;