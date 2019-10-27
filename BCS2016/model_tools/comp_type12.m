function X=comp_type12(pr,dum,Dsc,th)
% Dcm: dummy variable
%      1 if model predicts a crisis in t+1 | no crisis in t
dcm = (pr>th).*(dum==0); 
% The model predicts correctly the absence of crisis
S1  = (dcm==0).*(Dsc==0);
% Type I error: the model fails to predict a crisis
N1  = (dcm==0).*(Dsc==1);
p1  = sum(N1)/sum(Dsc==1);
% The model predicts correctly a crisis
S2  = (dcm==1).*(Dsc==1);
% Type II error: the model mistakenly predicts a crisis
N2  = (dcm==1).*(Dsc==0);
p2  = sum(N2)/sum(Dsc==0);
% Signal Noise ratio
Ns  = sum(S2)./sum(N2);
X   = [p1;p2;sum(dcm)];