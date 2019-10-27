function as=compute_steady_state(approx,apprf,data,param)
% as=compute_steady_state(approx,apprf,data,param)
%
% Input:    approx  : Structure for approximation coefficients of a(t+1)
%           apprf   : Structure for approximation coefficients of riskfree rate
%           data    : Structure that passes information for the decision
%                     rules (see main file)
%           param   : Parameters
%
% Output:   as      : Cell structure of nz elements, each containing the
%                     steady state associated with each level of the technology shock
%           Example : as{12}.y contains the steady state output level for the
%                     12th value of the technology shock.
%
%% Computes steady states
lb      = 0.5;
a0      = data.kss;
itmax   = 2000;
tol     = data.tol;
nz      = size(data.Gz,1);
as      = cell(nz,1);
stst    = zeros(nz,1);
for i=1:nz
    crit    = 1;
    iter    = 1;
    while and(crit>tol,iter<=itmax)
        a1  = comp_ap(a0,data.abar,approx,i);
        crit= abs(a1-a0);
        a0  = (1-lb)*a0+lb*a1;
        iter= iter+1;
    end
    
    if iter<itmax
        as{i}   = compute_dynamics(a0,i,data,approx,apprf,param,1);
        names   = fieldnames(as{i});
        stst(i) = 0;
        
        as{i}.rd	= comp_ap(as{i}.a1,data.abar,apprf,i);
        as{i}.dse   = param(2)./(as{i}.rd-param(2));
        as{i}.levb  = as{i}.dse+as{i}.phi*(1+as{i}.dse);
        as{i}.levl  = as{i}.dse;
        as{i}.leva  = as{i}.p^param(7)*as{i}.levl+(1-as{i}.p^param(7))*as{i}.levb;
        as{i}.sr    = 1-as{i}.c/(as{i}.x+as{i}.c);
    else
        stst(i) = 1;
        for j=1:length(names)
            as{i}.(names{j})   = NaN;
        end
    end
end

if any(stst)
    I=find(stst==1);
    fmt = '%d States admit no steady states:';
    fmt = strcat(fmt,repmat(' %d,',1,length(I)-1),' %d.\n');
    fprintf(fmt,length(I),I)
end