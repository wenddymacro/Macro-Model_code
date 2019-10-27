function [cc,coef]=cheb(xx,nn)
nx  = size(xx,1);
if nn>1;
    cc      = ones(nx,nn+1);
    coef    = zeros(nn+1,nn+1);
    cc(:,1) = 1;
    cc(:,2) = xx;
    coef(1:2,1:2)=eye(2);
    for i=3:nn+1;
        cc(:,i) = 2*xx.*cc(:,i-1)-cc(:,i-2);
        coef(i,1:i) = 2*[0 coef(i-1,1:i-1)]-[coef(i-2,1:i-2) 0 0];
    end
elseif nn==1;
    cc      = ones(nx,2);
    cc(:,2) = xx;
    coef    = eye(2);
else
    cc      = ones(nx,2);
    coef    = 1;
end
