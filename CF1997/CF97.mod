%=========================================================================%
% Carlstrom and Fuerst model (1997, AER)
% Ambrogio Cesa-Bianchi, November 2012
%=========================================================================%
% If you find any bugs when using this file or want to give me comments 
% and suggestions you can email me at ambrogio.cesabianchi@gmail.com

% I thank Kevin Salyer for spotting a bug in a previous version of this code


var K ke H He h q n i omegab ce cc w we Y r I Cc Ce C Bankrupcy Rb rpBANK rpENT lev rif PHI phi f g A;

varexo eA eN;

parameters alpha zeta beta delta gamma eta S nu mu p M rhoA StdeA StdeN;      
		alpha    = 0.36;  
		zeta     = 1-alpha-0.0001; 
		beta     = 0.99; 
		delta    = 0.02;  
		gamma    = 0.9474; 
		eta      = 0.1;   
		rhoA     = 0.99;
		StdeA    = 0.01;   
		StdeN    = 0.1;   
		S        = 0.207;  
		nu       = 2.52;   
		mu       = 0.25;   
		p        = pi;   
		M        = -.5*S^2;  
		rhoStdA  = 0.83;   
		StdeStdA = 0.19; 


model;
    q/cc = beta*(1/cc(+1))*(q(+1)*(1-delta)+r(+1));
    nu*cc = w;
    ke = i*f - ce/q;
    n = we + (ke(-1))*(q*(1-delta)+r) + StdeN*eN;
    K = (1-delta)*K(-1) + I*(1-mu*PHI);
    (1-eta)*cc + eta*ce + eta*i = Y;
    q =  (beta*gamma) * (r(+1) + q(+1)*(1-delta))  * ( (q(+1)*f(+1))/(1-q(+1)*g(+1)) );
    q = 1/(1 - mu*PHI - (mu*phi*f)/(1-PHI));
    i = (1/(1-q*g)) * n;
    Y = A*(K(-1)^alpha)*(H^zeta)*(He^(1-alpha-zeta));
    r = alpha*A*(K(-1)^(alpha-1))*(H^zeta)*(He^(1-alpha-zeta));
    w = zeta*A*(K(-1)^alpha)*(H^(zeta-1))*(He^(1-alpha-zeta));
    we = (1-alpha-zeta)*A*(K(-1)^alpha)*(H^zeta)*He^(-(alpha+zeta));
    H = (1-eta)*h;
    He = eta;
    Cc = (1-eta)*cc;       
    Ce = eta*ce;            
    C = (1-eta)*cc + eta*ce; 
    I = eta*i;              
    Bankrupcy = PHI;         
    Rb = (q*i*omegab)/(i-n); 
    rpBANK = Rb - 1;         
    lev = i/n;              
    rif = q*f*i/n;          
    rpENT = q*(1+r)-Rb;      
    PHI = normcdf((log(omegab)-M)/S);  
    phi = normpdf((log(omegab)-M)/S) / (omegab*S);                                    
    g   = normcdf((log(omegab)-M)/S - S) - PHI*mu + (1-PHI)*omegab; 
    f   = 1 - mu*PHI - g;    
    A    = (1-rhoA) + rhoA*A(-1) + StdeA*eA;
end;


initval;
    omegab  = 0.603892;
    PHI     = normcdf((log(omegab)-M)/S);
    phi     = normpdf((log(omegab)-M)/S) / (omegab*S);
    g       = normcdf((log(omegab)-M)/S - S) - PHI*mu + (1-PHI)*omegab;
    f       = 1-mu*PHI- g;
    q       = 1/(1-mu*0.01+(gamma-1)*f);
    r       = q*((1-beta*(1-delta))/beta);
    H       = .3;
    He      = eta;
    h       = .3/(1-eta);
    K       = (alpha/r)^(1/(1-alpha))*(He^(1-alpha-zeta))*(H^(zeta/(1-alpha)));
    Y       = (K^alpha)*(H^zeta)*(He^(1-alpha-zeta));
    i       = (delta/(eta*(1-mu*0.01)))*K;
    n       = (1-g*q)*i;
    ke      = (beta/q)*(eta*n-(1-alpha-zeta)*Y);
    ce      = q*(f*i-(ke/eta));
    cc      = (Y - eta*ce - eta*i)/(1-eta);
    A       = 1;
end;

steady;
check;

shocks;
    var eA = 1;
    var eN = 1; 
end;

stoch_simul(order=1,irf=24,nograph);


