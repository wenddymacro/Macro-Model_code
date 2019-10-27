% bgm_simp.mod
% Non-capital version of bilbiie, ghironi, melitz (2012)


var rhoo d w l le lc ne n v c z lambda k yc rk i;

varexo e;

predetermined_variables k;

parameters thetta betta fe chii psii phii deltaf r muu deltak alpha;

betta = .99;
deltaf = .025;
thetta = 3.8;
fe = 1;
chii = .924271;
psii = 4;
phii = .979;
deltak = .025;
alpha = .67;

r = betta^(-1)-1;

% Remains unchanged
muu = (thetta)/(thetta-1);
  
model;

% Changed - 1
rhoo = muu*lambda;

% Remains unchanged - 11
rhoo = n^(1/(thetta-1));

% Changed - 2
d = (1-1/muu)*(yc/n);

% Remains unchanged - 12
v = w*fe/z;

% Remains unchanged - 13
n = (1-deltaf)*(n(-1) + ne(-1));

% Remains unchanged - 14
% JR -- now changed. We need exogenous labor supply
% chii*l^(1/psii) = w/c;
l = 1;

% Remains unchanged - 15
v = betta*(1-deltaf)*(c/c(+1))*(v(+1)+d(+1));

% Changed - 5
% yc + ne*v = w*l + n*d + rk*k;

% Additional Equations

% Unchanged
ln(z) = phii*ln(z(-1)) + e;

% Unchanged - 10
le = ne*fe/z;

% Unchanged
lc = l - le;

% New - 3
k(+1) = (1-deltak)*k + i;

% New - 4
1=betta*(c(+1)/c)^(-1)*(rk(+1) + 1 - deltak);

% New - 6
yc = c + i;

% New - 7
w = (alpha/muu)*(yc/lc);

% New - 8
rk = ((1-alpha)/muu)*(yc/k);

% New - 9
yc = rhoo*z*lc^(alpha)*k^(1-alpha);

end;

%initval;
 steady_state_model;

z = 1;

l = 1;

n = (1-deltaf)*z*l/(fe*(alpha*(r+deltaf)/(muu-1) + deltaf));

rhoo = n^(1/(thetta-1));

lambda = rhoo/muu;

k = (z*(1-alpha)*rhoo/((r+deltak)*muu))^(1/alpha)*(l-n*(fe/z)*(deltaf/(1-deltaf)));

ne = (deltaf/(1-deltaf))*n;

le = ne*fe/z;

lc = l - le;

yc = rhoo*z*lc^(alpha)*k^(1-alpha);

i = deltak*k;

c = yc - i;

w = (alpha/muu)*(yc/lc);

rk = ((1-alpha)/muu)*(yc/k);

d = (1-1/muu)*(yc/n);

v = w*fe/z;


end;

steady;
check(qz_zero_threshold=1e-20);


shocks;
var e = .0072^2;

end;



stoch_simul (order=2,irf=0)  ;



