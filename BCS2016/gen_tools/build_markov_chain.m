function [Ga,Pa]= build_markov_chain(ra,sea,na,app)
%
% Builds a markov chain for an AR(1) process
% 
% [Ga,Pa]= build_markov_chain(ra,sea,na,type)
%
%  Inputs:  ra  : persistence of AR(1)
%           sea : std. dev of innovation
%           na  : Number of states
%           app : Approximation Method, 1 -> Rouwenhorst approximation
%                                       0 -> Tauchen Hussey
%
switch app
    case 1
        sa      = sea/sqrt(1-ra*ra);
        pa      = (1+ra)/2;         
        Pa      = rouwenhorst(na,pa);
        Ga      = exp(linspace(-sqrt(na-1)*sa,sqrt(na-1)*sa,na))';
    otherwise
        [xx,wx] = gauss_herm(na);       
        Ga      = exp(sqrt(2)*sea*xx);  
        x       = xx(:,ones(na,1));
        y       = x';
        w       = wx(:,ones(na,1))';
        %
        % computation
        %
        Pa      = (exp(y.*y-(y-ra*x).*(y-ra*x)).*w)./sqrt(pi);
        sx      = sum(Pa,2);
        Pa      = Pa./sx(:,ones(na,1));
end