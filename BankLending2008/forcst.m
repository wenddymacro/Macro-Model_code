function [yf,int_width,int_width2]=forcst(dr,y0,horizon,var_list)

% function [yf,int_width]=forecst(dr,y0,horizon,var_list)
%   computes mean forecast for a given value of the parameters
%   computes also confidence band for the forecast
%
% INPUTS:
%   dr:          structure containing decision rules
%   y0:          initial values
%   horizon:     nbr of periods to forecast
%   var_list:    list of variables (character matrix)
%
% OUTPUTS:
%   yf:          mean forecast
%   int_width:   distance between upper bound and
%                mean forecast
%
% SPECIAL REQUIREMENTS
%    none
%
% part of DYNARE, copyright Dynare Team (2003-2008)
% Gnu Public License.

   global M_  oo_ options_

   make_ex_;
   yf = simult_(y0,dr,zeros(horizon,M_.exo_nbr),1);
   nstatic = dr.nstatic;
   npred = dr.npred;
   nc = size(dr.ghx,2);
   endo_nbr = M_.endo_nbr;
   inv_order_var = dr.inv_order_var;
   [A,B] = kalman_transition_matrix(dr,nstatic+(1:npred),1:nc,dr.transition_auxiliary_variables);

   nvar = size(var_list,1);
   if nvar == 0
       nvar = M_.endo_nbr;
       ivar = [1:nvar];
   else
       ivar=zeros(nvar,1);
       for i=1:nvar
           i_tmp = strmatch(var_list(i,:),M_.endo_names,'exact');
           if isempty(i_tmp)
               disp(var_list(i,:));
               error (['One of the variable specified does not exist']) ;
           else
               ivar(i) = i_tmp;
           end
       end
   end

   ghx1 = dr.ghx(inv_order_var(ivar),:);
   ghu1 = dr.ghu(inv_order_var(ivar),:);

   sigma_u = B*M_.Sigma_e*B';
   sigma_u1 = ghu1*M_.Sigma_e*ghu1';
   sigma_y = 0;

   for i=1:horizon
       sigma_y1 = ghx1*sigma_y*ghx1'+sigma_u1;
       var_yf(i,:) = diag(sigma_y1)';
       if i == horizon
           break
       end
       sigma_u = A*sigma_u*A';
       sigma_y = sigma_y+sigma_u;
   end

   fact = qnorm((1-options_.conf_sig)/2,0,1);
   fact2 = qnorm((1-0.5)/2,0,1);

   int_width = zeros(horizon,M_.endo_nbr);
   int_width2 = zeros(horizon,M_.endo_nbr);
   for i=1:nvar
       int_width(:,i) = fact*sqrt(var_yf(:,i));
       int_width2(:,i) = fact2*sqrt(var_yf(:,i));
   end

   yf = yf(ivar,:);
   