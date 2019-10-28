---------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var x R pi g z xobs Robs pieobs Epi w;
varexo e_R e_g e_z e_v;

parameters Tau kappa rho_R rho_g rho_z psi1 psi2 rstar piestar ;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------



%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model(linear);

# beta = (1+(rstar/100))^(-1/4);
# tau = 1/Tau;

x = x(+1) - tau*(R - pi(+1)) + g;

pi = beta*pi(+1) + kappa*(x - z);

R = rho_R*R(-1) + (1-rho_R)*(psi1*pi + psi2*(x-z)) + e_R;
g = rho_g*g(-1) + e_g;
z = rho_z*z(-1) + e_z;

Epi=pi(+1);
w=1/(psi1+(1-beta)/kappa*psi2)*w(-1)+e_v-(pi-Epi(-1));


% Observation equations:
xobs = x;
pieobs = piestar + 4*pi;
Robs = rstar + piestar + 4*R;

end;



estimated_params;

psi1,1.2, , ,gamma_pdf,1.1,0.5;
psi2, , , ,gamma_pdf,0.25,0.15;
rho_R,  , , ,beta_pdf,0.5,0.2;
piestar,  , , ,gamma_pdf,4,2;
rstar,  , , ,gamma_pdf,2,1;
kappa, , , ,gamma_pdf,0.5,0.35;
Tau,  , , ,gamma_pdf,2,0.5;
rho_g,   , , , beta_pdf, 0.7, 0.1;
rho_z,  , , , beta_pdf, 0.7, 0.1;
stderr e_R,  , , , inv_gamma_pdf, 0.31, 0.16;
stderr e_g,  , , , inv_gamma_pdf, 0.38, 0.2;
stderr e_z,  , , , inv_gamma_pdf, 1, 0.52;
corr e_z, e_g, , -1, 1, normal_pdf, 0, 0.4;
corr e_R, e_g, , -1, 1, normal_pdf, 0, 0.4;
corr e_z, e_R, , -1, 1, normal_pdf, 0, 0.4;

stderr e_v,  , , , inv_gamma_pdf, 2, 0.5;
corr e_R, e_v, , -1, 1, normal_pdf, 0, 0.1;
corr e_g, e_v, , -1, 1, normal_pdf, 0, 0.1;
corr e_z, e_v, , -1, 1, normal_pdf, -0.7, 0.1;

end;


varobs xobs Robs pieobs;



estimation(datafile=nicoloindeter, xls_range=A1:C501,mode_compute=0,mode_file=augmentedrep_mode,mh_replic=250000,mh_drop=0.5, mh_jscale=0.4,graph_format=(pdf,fig),mode_check);
