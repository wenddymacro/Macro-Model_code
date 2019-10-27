% bgm_simp.mod
% Non-capital version of bilbiie, ghironi, melitz (2012)


var vc vy vn vnc vny Y YI rhoo d w l le lc ne n v c z ln_ca ln_c ln_ne ln_n ln_d ln_da ln_v ln_vva ln_l ln_le ln_lc ln_w ln_wa ln_rhoo ln_Y ln_YA ;

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

le = ne*fe/z;

lc = l - le;

Y = c + ne*v;
YI = w*l + n*d;

% Empirically relevant, log

ln_c = log(c);
ln_ca = log(c/rhoo);

ln_d = log(d);
ln_da = log(d/rhoo);

ln_w = log(w);
ln_wa = log(w/rhoo);

ln_v = log(v);
ln_vva = log(v/rhoo);

ln_Y = log(Y);
ln_YA = log(Y/rhoo);

ln_l = log(l);

ln_le = log(le);
ln_lc = log(lc);

ln_ne = log(ne);
ln_n = log(n);

ln_rhoo = log(rhoo);

vc = v/c;
vy = v/Y;

vn = v*n;

vnc = vn/c;
vny = vn/Y;

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

Y = c + ne*v;
YI = w*l + n*d;

le = ne*fe/z;

lc = l - le;

ln_ca = log(c/rhoo);
ln_c = log(c);

ln_d = log(d);
ln_da = log(d/rhoo);

ln_v = log(v);
ln_vva = log(v/rhoo);

ln_w = log(w);
ln_wa = log(w/rhoo);

ln_Y = log(Y);
ln_YA = log(Y/rhoo);

ln_l = log(l);
ln_le = log(le);
ln_lc = log(lc);

ln_ne = log(ne);
ln_n = log(n);

ln_rhoo = log(rhoo);

vc = v/c;
vy = v/Y;
vn = v*n;

vnc = vn/c;
vny = vn/Y;

end;

steady;
check;


shocks;
var e = .0072^2;

end;



stoch_simul (order=2,irf=25) ln_n ln_ne  ;



