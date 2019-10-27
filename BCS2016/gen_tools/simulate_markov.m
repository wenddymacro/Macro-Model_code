function [chain,st] = simulate_markov(x,P,T,varargin)
%  x     = the quantity corresponding to each state, typical element x(i)
%  P     = Markov transition matrix, typical element p(i,j) i,j=1,...n
%  T     = number of periods to simulate
%  varargin{1}=start (to impose a particular initial state)
%  varargin{2}= seed
%  chain = sequence of realizations from the simulation
%  Modification of progam by L&S.
% if nargin>4;
%     rand('state',varargin{2})
% end
x   = x(:);
n   = size(x,1);    % what is the size of the state vector?
E   = rand(T,1);    % T-vector of draws from independent uniform [0,1]

cumsumP = P*triu(ones(size(P)));
%% creates a matrix whose rows are the cumulative sums of
%% the rows of P

%%%%% SET INITIAL STATE USING STATIONARY DISTRIBUTION
if and(nargin>3,~isempty(varargin{1}))
    nn  = length(varargin{1});
    s   = zeros(n,nn);
    for i=1:nn
        s(varargin{1}(i),i)=1;
    end
else
    nn      = 1;
    [Q,D]   = eig(P');
    pi0     = Q(:,(abs(diag(D)-1)<1e-12))';
    pi0     = pi0/sum(pi0);
    E0      = rand(1,1);
    ppi0    = [0,cumsum(pi0)];
    s0      = ((E0<=ppi0(2:n+1)).*(E0>ppi0(1:n)))';
    s       = s0;
end


%%%%% ITERATE ON THE CHAIN
state   = zeros(n,T);
st      = zeros(T,1);
state(:,1:nn)   = s;
[i,j]           = find(s);
st(1:nn)        = i;
s               = s(:,nn);
for t=nn:T,    
    state(:,t) = s;
    st(t)      = find(s);
    ppi        = [0,s'*cumsumP];
    s          = (+and(E(t)<=ppi(2:n+1),E(t)>ppi(1:n)))';
end
chain   = (x'*state)';


