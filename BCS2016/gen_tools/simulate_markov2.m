function [chain1,chain2] = simulate_markov2(x,P,T,varargin)
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
%%%%% SET INITIAL STATE USING STATIONARY DISTRIBUTION
if nargin>3
    init= varargin{1};
    nn  = size(init,2);
    s1  = zeros(n,nn);
    s2  = zeros(n,nn);
    for i=1:nn
        s1(init(1,i),i)=1;
        s2(init(2,i),i)=1;
    end
else
    nn      = 1;
    [Q,D]   = eig(P');
    pi0     = Q(:,(abs(diag(D)-1)<1e-12))';
    pi0     = pi0/sum(pi0);
    E0      = rand(1,1);
    ppi0    = [0,cumsum(pi0)];
%     s0      = ((E0<=ppi0(2:n+1)).*(E0>ppi0(1:n)))';
    s0      = ((E0<=ppi0(2:n)).*(E0>ppi0(1:n-1)))';
    s1      = [0;s0];
    s2      = [s0;0];
end

%%%%% ITERATE ON THE CHAIN
state1          = zeros(n,T);
state2          = zeros(n,T);
st1             = zeros(T,1);
st2             = zeros(T,1);
state1(:,1:nn)  = s1;
state2(:,1:nn)  = s2;
[i1,~]          = find(s1);
[i2,~]          = find(s2);
st1(1:nn)       = i1;
st2(1:nn)       = i2;
s1              = s1(:,nn);
s2              = s2(:,nn);
for t=nn:T,
    state1(:,t) = s1;
    state2(:,t) = s2;
    st1(t)      = find(s1);
    st2(t)      = find(s2);
    ppi1        = [0,s1'*cumsumP];
    ppi2        = [0,s2'*cumsumP];
    s1          = (+and(E(t)<=ppi1(2:n+1),E(t)>ppi1(1:n)))';
    s2          = (+and(E(t)<=ppi2(2:n+1),E(t)>ppi2(1:n)))';
end
chain1  = (x'*state1)';
chain2  = (x'*state2)';


