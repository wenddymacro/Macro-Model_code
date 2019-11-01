% Generate probability distributions of A0, A+, or impulse responses
%       See WZ's ``A Gibbs sampler for structural VARs.''
%
% August 2000, Tao Zha


msstart2    % start the program in which everyhting is initialized through mpstart.m
if ~indxEstima
   warning('You must set indxEstima=1 in msstart to run this program')
	disp('Press ctrl-c to abort now')
	pause
end

imfgbs_draws = zeros(ndraws2,imstp*nvar^2);  % Storing draws of impulse responses.
FM_g_gbs_draws = zeros(ndraws2,imstp);  %Storing draws of finanical multiplier of government spending.
FM_t_gbs_draws = zeros(ndraws2,imstp);  %Storing draws of finanical multiplier of government taxes.

%if indxModStr  %nCsk    % manually change case by case.  Conditional on structural shocks
%           % So far this option contains a fatal error because eplhate should change
%           %   with each draw of A0 and A+.  09/15/00.
%   %==== Get the last row of X -- phil -- for conditional forecasting
%   Bfcrt = [2005 3];   % beginning of the forecast dates.
%   Efcrt = [2007 3];  %[yrEnd qmEnd];  % end of the forecast dates.  Note must be earlier than yrEnd and qmEnd.
%   Estre = fn_dataext(Bfcrt,Efcrt,eplhate);  % backed out structural shocks with dates
%   tmp = fn_dataext(Bfcrt,Bfcrt,phie);
%   phil = tmp(:,[3:end]);
%   %
%   Estr1 = Estrcon(:,3:end);
%%   Estr1 = Estre(:,3:end);      % used for conditional forecasting later
%%   Estr1(:,[1 3:7])=zeros(size(Estre(:,1)),6);   % only MS shocks
%   %Estr1(:,1:8)=zeros(size(Estre(:,1)),8);
%end

%-------- Preallocating space to store distributions -----------------
%--- Note: except *mean and *hat for imf and yr, all matrices are vectorized for imf and yr
%
bmean = zeros(sum(n0),1);  % free parameters b -- stacked vector of free parameters in A0.
gmean = zeros(sum(np),1);  % free parameters g -- stacked vector of free parameters in A+.
imfmean = zeros(imstp,nvar^2);
                       % nvar responses to 1st shock, nvar responses to 2nd shock, etc.
yrgumean  = zeros(nfyr,nvar);  % unconditional forecasts: annual rates that work for
           % both calendar and any annual years, depending on Byrqm and Eyrqm
           % with reshape(ygmean,nfyr,nvar) if needed to be in the 2-d matrix form.
qmyrgumean = zeros(nfqm,nvar);  % unconditional forecasts: period-to-last-period annual growth rates
           % with reshape(qmygmean,nfqm,nvar) if needed to be in the 2-d matrix form.
%
bcnt = zeros(ninv+2,sum(n0));  % count (dist) of free parameters b -- stacked vector of free parameters in A0.
gcnt = zeros(ninv+2,sum(np));  % count (dist) of free parameters g -- stacked vector of free parameters in A+.
imfcnt = zeros(ninv+2,imstp*nvar^2);  % impulse responses
yrgucnt  = zeros(ninv+2,nfyr*nvar);  % unconditional annual rates, works for both calendar and any annual years, depending on Byrqm and Eyrqm
qmyrgucnt = zeros(ninv+2,nfqm*nvar);  % unconditional period-to-last-period annual growth rates


%=================================================
% Generate error bands by first drawing from the
%  t-distribution for each column equation
%=================================================
%
%*** The following may not be used elsewhere
%  bhat = cell(nvar,1);  % ML b: cell i represents free parameters b(i) in the ith equation
%  n0cum = [0;cumsum(n0)];
%  for kj = 1:nvar
%     bhat{kj} = xhat(n0cum(kj)+1:n0cum(kj+1));
%  end
% bgbs:  new draw of b -- a stacked sum(n0)-by-1 vector of free parameters in A0

%~~~~~~~~~~~~ The following is for the cross-A0-A+ restrictions ~~~~~~~~~~~~
if indxC0Pres
   ixeqv=zeros(size(ixmC0Pres,1),1);  % v: vector; eq: equation
   Lit = cell(size(ixmC0Pres,1),1);  % transposed restriction matrix Li, which is
                                    % V_j(i,:) in f_j(i) = V_j(i,:)*g_j
   LtH = cell(size(ixmC0Pres,1),1);   % Li'*H
   LV = cell(size(ixmC0Pres,1),1);   % Li*inv(LtH*Li)
   stdR = cell(size(ixmC0Pres,1),1);    % lower triagular chol of restricted covariance
                  % see Zha's forecast (1), p.17
   for ki = 1:size(ixmC0Pres,1)   % loop through the number of equations in which
               % cross-A0-A+ restrictions occur. See St. Louis Note p.5.
      ixeqv(ki) = ixmC0Pres{ki}(1,1);  % should be all in the set of 1:nvar;
         %  indecies for these equations, e.g., [1 3]' means 1st and 3rd equations.
      Lit{ki} = Vi{ixeqv(ki)}(ixmC0Pres{ki}(:,2),:);  % transposed restriction matrix Li
               % V_j(i,:) in f_j(i) = V_j(i,:)*g_j
      LtH{ki} = Lit{ki}/Hpinv{ixeqv(ki)};   % Li'*H
      LV{ki} = Lit{ki}'/(LtH{ki}*Lit{ki}');   % Li*inv(LtH*Li)
      covR = Hpinv{ixeqv(ki)}\( eye(size(Hpinv{ixeqv(ki)}))-LV{ki}*LtH{ki} );
               % Restricted covariance Hi- Hi*Li*inv(Li'*Hi*Li)*Li'*Hi
               % May be singular but non-negative, symmetric.
      [u1 d1] = eig(covR);
      %[u1 d1 v1] = svd(Ome);  % too slow
      stdR{ki} = u1*diag(sqrt(diag(abs(d1))));    % lower triagular chol of restricted covariance
                     % see Zha's forecast (1), p.17
   end
end


if IndxParR    % with paramaters random, starting with Gibbs draws of A0 and inv(A0)

   %(1)---------- Preparations outside the Gibbs loop ------------
   %**** For A0's
   n0cum = [0;cumsum(n0)];
   H0sqr_inv = zeros(n0cum(end));  % inv(sqrt(H0)) where H) has NOT divided by T.
   bhat=xhat;   % ML estimate of b -- a stacked sum(n0)-by-1 vector of free parameters in A0
   %
   A0gbs = A0hat;   % A0's from the Gibbs sampling
   cTinv = cell(nvar,1);   % in each cell, inv(T_i) for T_iT_i'=H0_i in the proof of Theorem 1 in the WZ paper
   UT = cell(nvar,1);  % in each cell, U_i*T_i in the proof of Theorem 1
   for k=1:nvar
      cTinv{k} = chol(H0inv{k}/fss);   % upper triangular Choleski for H0inv but lower triangular Choleski to H0inv
      UT{k} = Ui{k}/cTinv{k};    % n-by-qi: U_i*T_i in the proof of Theorem 1
      %
      H0sqr_inv(n0cum(k)+1:n0cum(k+1),n0cum(k)+1:n0cum(k+1)) = cTinv{k}*sqrt(fss);
   end

   %**** For A+'s
   VHphalf = cell(nvar,1);  % in each cell, V_i*sqrt(Hp_i).
   PU = cell(nvar,1);  % in each cell, P_i*U_i
   VPU = cell(nvar,1);  % in each cell, V_i*P_i*U_i
   for ki=1:nvar
      VHphalf{ki} = Vi{ki}/chol(Hpinv{ki}); % where chol(Hpinv_i)*chol(Hpinv_i)'=Hpinv_i.
      PU{ki} = Pmat{ki}*Ui{ki}';
      VPU{ki} = Vi{ki}*PU{ki};
   end

   %*** For forecasts (conditional or unconditional) without including dates in the first 2 columns
   if (E2yrqm(end,2)~=q_m)
      temp = fn_dataext([fdates(1,1) 0],[fdates(end-E2yrqm(end,2),1) 0],yafyrghate);
   else
      temp = fn_dataext([fdates(1,1) 0],[fdates(end-1,1) 0],yafyrghate);
   end
   yrguhat = temp(:,3:end);   % nfry-by-nvar ML unconditional forecasts (annual growth)
   temp = fn_dataext(fdates(1,:),fdates(end,:),yafqmyghate);
   qmyrguhat = temp(:,3:end);
        % nfqm-by-nvar ML unconditional forecasts (annualized growth -- to q or m in the last year)


   %(2)------------------- Gibbs sampling ---------------------
   if (1)
      %*** Initial draw of A0 at ML
      A0gbs = A0hat;   % from "load ..."

      Fgbs = Fhat;  % initialize A+
   else
      %*** Initial drawe of A0 other than ML
      bgbs = H0sqr_inv\randn(sum(n0),1);  %H_sr*randn(nfp,1);
      csq=randn(tdf,1);
      csq=sum(csq .* csq);
      bgbs = bhat+bgbs/sqrt(csq/tdf);

      %** Normalization by the choice of IndxNmlr
      A0gbs = fn_tran_b2a(bgbs,Ui,nvar,n0);
      A0gbs = nmlzvar(A0gbs,A0hat,[],IndxNmlr,0,[]);

      Fgbs = randn(ncoef,nvar);  % initialize A+
   end


   %*** First set of draws of A0 for two distinct purposes:
   %           (1) to be discarded in the Gibbs sampling
   %           (2) to get the standard deviations of variables of interest
   b1=zeros(sum(n0),1); b2=b1; bstd=NaN; bRange=NaN; bhbin=NaN;
   g1=zeros(sum(np),1); g2=g1; gstd=NaN; gRange=NaN; ghbin=NaN;
   imf1=zeros(imstp*nvar^2,1); imf2=imf1; imfstd=NaN; imfRange=NaN; imfhbin=NaN;
   yrgu1=zeros(nfyr*nvar,1); yrgu2=yrgu1; yrgustd=NaN; yrguRange=NaN; yrguhbin=NaN;
   qmyrgu1=zeros(nfqm*nvar,1); qmyrgu2=qmyrgu1; qmyrgustd=NaN; qmyrguRange=NaN; qmyrguhbin=NaN;
   for draws = 1:ndraws1
      %*** A Gibbs draw of A0
      A0gbs = fn_gibbsrvar(A0gbs,UT,nvar,fss,n0,Indxcol);
      A0gbs = nmlzvar(A0gbs,A0hat,[],IndxNmlr,0,[]);  % normalization
      %
      if Aband   % get std for b or A0
         bgbs = fn_tran_a2b(A0gbs,Ui,nvar,n0);   % converted to b from A0
         b1 = b1+bgbs;     % 1st moment
         b2 = b2+bgbs.^2;  % 2nd moment
      end

      if IndxAp   % A+, impulse responses, and forecasts are all in this if loop.
         %*** A normal draw of A+
         for kj=Indxcol
            if indxC0Pres   % under cross-A0-A+ restrictions
               if isempty(find(kj==ixeqv))   % no cross-A0-A+ restrictions for the kj_th equation
                  Fgbs(:,kj) = VPU{kj}*A0gbs(:,kj) + VHphalf{kj}*randn(np(kj),1);   % A+
               else
                  ci = ixmC0Pres{kj}(:,4) .* A0gbs(ixmC0Pres{kj}(:,3),ixeqv(kj));
                           % s * a_j(h) in the restriction f_j(i) = s * a_j(h).
                  gbar = PU{kj}*A0gbs(:,kj);
                  fgbsmean = VPU{kj}*A0gbs(:,kj) + Vi{kj}*(Hpinv{kj}\( LV{kj}*(ci-Lit{kj}*gbar) ));
                  Fgbs(:,kj) = fgbsmean + Vi{kj}*stdR{kj}*randn(np(kj),1);   % A+
                                     % see Zha's St. Louis Note, p.5
               end
            else     % no cross-A0-A+ restrictions
               Fgbs(:,kj) = VPU{kj}*A0gbs(:,kj) + VHphalf{kj}*randn(np(kj),1);   % A+
            end
         end
         %
         if Apband  % get std for g or A+
            ggbs = fn_tran_f2g(Fgbs,Vi,nvar,np); % a stacked vector of all free parameters in A+.
            g1 = g1+ggbs;    % 1st moment
            g2 = g2+ggbs.^2; % 2nd moment
         end

         if IndxImf   % draws of impulse responses
            swish = inv(A0gbs);
            Bgbs = Fgbs*swish;   % ncoef-by-nvar reduced form lagged parameters.
            imfgbs = fn_impulse(Bgbs,swish,[nvar lags imstp]);
                        % imstp-by-nvar^2 in the form that is congenial to RATS and
                        % in the order of nvar responses to 1st shock; nvar responses to 2nd shock, etc.
            if Imfband    % error bands or distributions of impulse responses
               imf1 = imf1+imfgbs(:);    % 1st moment
               imf2 = imf2+imfgbs(:).^2;    % 2nd moment
            end
         end
         %
         if IndxFore   % draws of out-of-sample forecasts
            %
            if ~IndxImf   % Impulse responses compute Bgbs.  Without them, Bgbs has to be computed in the loop.
               Bgbs = Fgbs/A0gbs;   % ncoef-by-nvar reduced form lagged parameters.
            end

            if nCsk     % conditional structural shocks
                        % So far this option contains a fatal error because eplhate should change
                        %   with each draw of A0 and A+.  09/15/00.
               nn = [nvar lags size(Estre,1)];
               if nexo<2
                  yforegbs = fn_forecastfixe(Bhat,A0hat,phil,nn,Estr1);  % *-by-nvar, in log
               else
                  Xfexoe = fn_dataext(Bfcrt,Efcrt,xdatae(:,[1:2 2+nvar+1:2+nvar+nexo-1]));
                  yforegbs = fn_forecastfixe(Bhat,A0hat,phil,nn,Estr,nexo,Xfexoe(:,3:end));    % *-by-nvar, in log
               end
            else     % without fixed structural shocks -- i.e., unconditional
               if nexo<2
                  yforegbs = fn_forecastsim(Bgbs,A0gbs,phil,[nvar lags nfqm]);
                                 % unconditional with simulated shocks; nfqm-by-nvar, in log
               else
                  yforegbs = fn_forecastsim(Bgbs,A0gbs,phil,[nvar lags nfqm],nexo,Xfexoe(:,3:end));
                                 % unconditional with simulated shocks; nfqm-by-nvar, in log
               end
            end
            yforegbse = [fdates yforegbs];

            yafgbse = [yact1e; yforegbse];  % actual and unconditional forecast.  Need actual to compute annual growth rates
            %===== Converted to mg (or qg) and yg
            [yafyrggbse,jnk,yafqmyggbse] = fn_datana(yafgbse,q_m,vlistlog(1:nlogeno),vlistper(1:npereno));
                              % actual and unconditional forecast growth rates

            if (E2yrqm(end,2)~=q_m)
               yforeyrggbse = fn_dataext([fdates(1,1) 0],[fdates(end-E2yrqm(end,2),1) 0],yafyrggbse);
            else
               yforeyrggbse = fn_dataext([fdates(1,1) 0],[fdates(end,1) 0],yafyrggbse);
            end
            yforeqmyggbse = fn_dataext(fdates(1,:),fdates(end,:),yafqmyggbse);
                                      % unconditional forecast growth rates with the dates
            if Foreband
               temp = yforeyrggbse(:,3:end);
               yrgu1 = yrgu1+temp(:);     % 1st moment  annual growth rate in unconditional forecast
               yrgu2 = yrgu2+temp(:).^2;  % 2nd moment  annual growth rate in unconditional forecast
               temp = yforeqmyggbse(:,3:end);
               qmyrgu1 = qmyrgu1+temp(:);  % 1st moment annualized (to q or m in the last year) growth rate in unconditional forecast
               qmyrgu2 = qmyrgu2+temp(:).^2;  % 2nd moment annualized (to q or m in the last year) growth rate in unconditional forecast
            end
         end
         %
      end
   end

   if Aband
      b1 = b1/ndraws1;
      b2 = b2/ndraws1;
      bstd = abs(sqrt(b2 - b1.^2));   % standard deviation.  Abs is used to gauard agaist almost zero but imaginary numbers.
      bstd(find(bstd<sqrt(eps)))=1e-04;   % set to a very small number for zero standard deviations.
      bRange=zeros(sum(n0),2);  % 2: min and max
      bRange(:,1) = bhat-nstd*abs(bstd);
      bRange(:,2) = bhat+nstd*abs(bstd);
      bhbin = (bRange(:,2) - bRange(:,1)) ./ ninv;
   end
   %
   if Apband
      g1 = g1/ndraws1;
      g2 = g2/ndraws1;
      gstd = abs(sqrt(g2 - g1.^2));   % standard deviation.  Abs is used to gauard agaist almost zero but imaginary numbers.
      gstd(find(gstd<sqrt(eps)))=1e-04;   % set to a very small number for zero standard deviations.
      gRange=zeros(sum(np),2);  % 2: min and max
      gRange(:,1) = ghat-nstd*abs(gstd);
      gRange(:,2) = ghat+nstd*abs(gstd);
      ghbin = (gRange(:,2) - gRange(:,1)) ./ ninv;

   end
   %
   if IndxAp
      if IndxImf
         if Imfband
            imf1 = imf1/ndraws1;
            imf2 = imf2/ndraws1;
            imfstd = abs(sqrt(imf2 - imf1.^2));   % standard deviation.  Abs is used to gauard agaist almost zero but imaginary numbers.
            imfstd(find(imfstd<sqrt(eps)))=1e-04;   % set to a very small number for zero standard deviations.
            imfRange=zeros(imstp*nvar^2,2);  % 2: min and max
            imfRange(:,1) = imfhat(:)-nstd*abs(imfstd);
            imfRange(:,2) = imfhat(:)+nstd*abs(imfstd);
            imfhbin = (imfRange(:,2) - imfRange(:,1)) ./ ninv;
         end
      end
      %
      if IndxFore
         if Foreband
            yrgu1 = yrgu1/ndraws1;
            yrgu2 = yrgu2/ndraws1;
            yrgustd = abs(sqrt(yrgu2 - yrgu1.^2));   % standard deviation.  Abs is used to gauard agaist almost zero but imaginary numbers.
            yrgustd(find(yrgustd<sqrt(eps)))=1e-04;   % set to a very small number for zero standard deviations.
            yrguRange=zeros(nfyr*nvar,2);  % 2: min and max
            yrguRange(:,1) = yrguhat(:)-nstd*abs(yrgustd);
            yrguRange(:,2) = yrguhat(:)+nstd*abs(yrgustd);
            yrguhbin = (yrguRange(:,2) - yrguRange(:,1)) ./ ninv;
            %
            qmyrgu1 = qmyrgu1/ndraws1;
            qmyrgu2 = qmyrgu2/ndraws1;
            qmyrgustd = abs(sqrt(qmyrgu2 - qmyrgu1.^2));   % standard deviation.  Abs is used to gauard agaist almost zero but imaginary numbers.
            qmyrgustd(find(qmyrgustd<sqrt(eps)))=1e-04;   % set to a very small number for zero standard deviations.
            qmyrguRange=zeros(nfqm*nvar,2);  % 2: min and max
            qmyrguRange(:,1) = qmyrguhat(:)-nstd*abs(qmyrgustd);
            qmyrguRange(:,2) = qmyrguhat(:)+nstd*abs(qmyrgustd);
            qmyrguhbin = (qmyrguRange(:,2) - qmyrguRange(:,1)) ./ ninv;
         end
      end
   end

   %=================================
   %*** Second set of draws of A0 to be kept
   %=================================
   tic
   draws = 0;
   indx_failgbs = 0;
   while (draws<=ndraws2)
      draws = draws + 1;
      if ~mod(draws,nbuffer)
         draws
         %  fwriteid = fopen('outA0.bin','a');
         %  count = fwrite(fwriteid,A0hatw,'double');
         %  status = fclose('all');
      end

      %*** A Gibbs draw of A0
      A0gbs = fn_gibbsrvar(A0gbs,UT,nvar,fss,n0,Indxcol);
      A0gbs = nmlzvar(A0gbs,A0hat,[],IndxNmlr,0,[]);  % normalization
      %
      if Aband     % get distributions for b or A0
         bgbs = fn_tran_a2b(A0gbs,Ui,nvar,n0);   % converted to b from A0
         bmean = bmean + bgbs;   % mean
         bcnt = fn_empdfsort(bcnt,bgbs,bRange(:,1),bhbin,ninv);
      end

      if IndxAp   % A+, impulse responses, and forecasts are all in this if loop.
         %*** A normal draw of A+
         for kj=Indxcol
            if indxC0Pres   % under cross-A0-A+ restrictions
               if isempty(find(kj==ixeqv))   % no cross-A0-A+ restrictions for the kj_th equation
                  Fgbs(:,kj) = VPU{kj}*A0gbs(:,kj) + VHphalf{kj}*randn(np(kj),1);   % A+
               else
                  ci = ixmC0Pres{kj}(:,4) .* A0gbs(ixmC0Pres{kj}(:,3),ixeqv(kj));
                           % s * a_j(h) in the restriction f_j(i) = s * a_j(h).
                  gbar = PU{kj}*A0gbs(:,kj);
                  fgbsmean = VPU{kj}*A0gbs(:,kj) + Vi{kj}*(Hpinv{kj}\( LV{kj}*(ci-Lit{kj}*gbar) ));
                  Fgbs(:,kj) = fgbsmean + Vi{kj}*stdR{kj}*randn(np(kj),1);   % A+
                                     % see Zha's St. Louis Note, p.5
               end
            else     % no cross-A0-A+ restrictions
               Fgbs(:,kj) = VPU{kj}*A0gbs(:,kj) + VHphalf{kj}*randn(np(kj),1);   % A+
            end
         end
         %
         if Apband     % get distributions for g or A+
            ggbs = fn_tran_f2g(Fgbs,Vi,nvar,np);   % % convert A to a stacked vector of all free parameters in A+.
            gmean = gmean + ggbs;   % mean
            gcnt = fn_empdfsort(gcnt,ggbs,gRange(:,1),ghbin,ninv);  % count or distribution
         end

         if IndxImf   % draws of impulse responses
            %--- Original Gibbs draw.
            A0gbsinv = inv(A0gbs);
            Bgbs = Fgbs*A0gbsinv;   % ncoef-by-nvar reduced form lagged parameters.
            %****** Sign restrictions ******
            if (0)
               [Q_SRgbs, indx_failgbs] = ftd_SRestrictRWZalg_FP_Model2(A0gbsinv,Bgbs,nvar,lags,npers_signs,ndraws_terminate);
               if (~indx_failgbs)
                  A0gbs = A0gbs* Q_SRgbs';
                  A0gbsinv = Q_SRgbs*A0gbsinv;
                  Fgbs = Bgbs*A0gbs;
               else
                  warning('Failed for sign restrictions')
                  draws = draws - 1; %Resetting the draws (very important).
               end
            end

            if (~indx_failgbs)   %Success for sign restrictions.
               swish = A0gbsinv;
               %*** No need for this. Bgbs = Fgbs*swish;   % ncoef-by-nvar reduced form lagged parameters.
               imfgbs = fn_impulse(Bgbs,swish,[nvar lags imstp]);
                           % imstp-by-nvar^2 in the form that is congenial to RATS and
                           % in the order of nvar responses to 1st shock; nvar responses to 2nd shock, etc.
               imf3gbs = reshape(imfgbs,imstp,nvar,nvar);
                           % imf3: row--steps, column--nvar responses, 3rd dimension--nvar shocks

               %=== Storing the draws of impulse responses.
               imfgbs_draws(draws,:) = imfgbs(:)';
               %---
               if Imfband    % error bands or distributions of impulse responses
                  imfmean = imfmean + imfgbs;  % mean
                  imfcnt = fn_empdfsort(imfcnt,imfgbs,imfRange(:,1),imfhbin,ninv);  % count or distribution
               end
            end
         end


         if (~indx_failgbs)      %Success for sign restrictions.
            if IndxFore   % draws of out-of-sample forecasts
               %
               if ~IndxImf   % Impulse responses compute Bgbs.  Without them, Bgbs has to be computed in the loop.
                  Bgbs = Fgbs/A0gbs;   % ncoef-by-nvar reduced form lagged parameters.
               end
               if nexo<2
                  yforegbs = fn_forecastsim(Bgbs,A0gbs,phil,[nvar lags nfqm]);
                                 % unconditional with simulated shocks; nfqm-by-nvar, in log
               else
                  yforegbs = fn_forecastsim(Bgbs,A0gbs,phil,[nvar lags nfqm],nexo,Xfexoe(:,3:end));
                                 % unconditional with simulated shocks; nfqm-by-nvar, in log
               end
               yforegbse = [fdates yforegbs];

               yafgbse = [yact1e; yforegbse];  % actual and unconditional forecast.  Need actual to compute annual growth rates
               %===== Converted to mg (or qg) and yg
               [yafyrggbse,jnk,yafqmyggbse,jnk2,yafyqgbse, jnk3, yaflevelgbse] = fn_datana2(yafgbse,q_m,vlistlog(1:nlogeno),vlistper(1:npereno));
                                 % actual and unconditional forecast growth rates
               yaflevelgbse_draws(:,:,draws) = yaflevelgbse;
               yafyqgbse_draws(:,:,draws) = yafyqgbse;


               if (E2yrqm(end,2)~=q_m)
                  yforeyrggbse = fn_dataext([fdates(1,1) 0],[fdates(end-E2yrqm(end,2),1) 0],yafyrggbse);
               else
                  yforeyrggbse = fn_dataext([fdates(1,1) 0],[fdates(end,1) 0],yafyrggbse);
               end
               yforeqmyggbse = fn_dataext(fdates(1,:),fdates(end,:),yafqmyggbse);
                                         % unconditional forecast growth rates with the dates
               if Foreband
                  yrgugbs = yforeyrggbse(:,3:end);
                  yrgumean = yrgumean+yrgugbs;     % mean.   nfyr-by-nvar
                  yrgucnt = fn_empdfsort(yrgucnt,yrgugbs,yrguRange(:,1),yrguhbin,ninv);
                                % count or distribution.   ninv+2-by-nfyr*nvar
                  %
                  qmyrgugbs = yforeqmyggbse(:,3:end);
                  qmyrgumean = qmyrgumean+qmyrgugbs;     % mean.  nfqm-by-nvar
                  qmyrgucnt = fn_empdfsort(qmyrgucnt,qmyrgugbs,qmyrguRange(:,1),qmyrguhbin,ninv);
                                % count or distribution.  ninv+2-by-nfqm*nvar
               end
            end
         end
      end
   end
   timeminutes=toc/60

   if Aband
      bmean = bmean/ndraws2;
   end
   %
   if Apband
      gmean = gmean/ndraws2;
   end
   %
   if IndxAp
      if IndxImf
         if Imfband
            imfmean = imfmean/ndraws2;
         end
      end
      %
      if IndxFore
         if Foreband
            yrgumean=yrgumean/ndraws2;
            qmyrgumean=qmyrgumean/ndraws2;
         end
      end
   end
else     % without paramaters random and conditional on, say, ML estimates.  Simulating error bands for forecasts.
         % In this section, I haven't got time to do forecasts wtih cross-A0-A+ restrictions.  9/18/00.

   %(1)---------- Preparations outside the simulation loop ------------
   %**** For A0's
   bhat=xhat;   % This is never used, but ``save outdistcnt ..'' calls this variable

   %*** For forecasts (conditional or unconditional) without including dates in the first 2 columns
   temp = fn_dataext([fdates(1,1) 0],[fdates(end,1) 0],yafyrghate);
   yrguhat = temp(:,3:end);   % nfry-by-nvar ML unconditional forecasts (annual growth)
   temp = fn_dataext(fdates(1,:),fdates(end,:),yafqmyghate);
   qmyrguhat = temp(:,3:end);
        % nfqm-by-nvar ML unconditional forecasts (annualized growth -- to q or m in the last year)


   %(2)------------------- Simulation loop ---------------------
   %*** The following 3 lines are of no use, but we have them initialized so that they
   %***   can be saved in ``save outdistcnt ..'' without causing the program to quit.
   b1=zeros(sum(n0),1); b2=b1; bstd=NaN; bRange=NaN; bhbin=NaN;
   g1=zeros(sum(np),1); g2=g1; gstd=NaN; gRange=NaN; ghbin=NaN;
   imf1=zeros(imstp*nvar^2,1); imf2=imf1; imfstd=NaN; imfRange=NaN; imfhbin=NaN;

   yrgu1=zeros(nfyr*nvar,1); yrgu2=yrgu1; yrgustd=NaN; yrguRange=NaN; yrguhbin=NaN;
   qmyrgu1=zeros(nfqm*nvar,1); qmyrgu2=qmyrgu1; qmyrgustd=NaN; qmyrguRange=NaN; qmyrguhbin=NaN;
   %*** First set of draws of forecasts for the purpose of
   %           of getting the standard deviations of variables of interest
   for draws = 1:ndraws1
      if IndxFore   % draws of out-of-sample forecasts
         if nconstr %WZ-DSL conditional forecasts.
            [ycforegbs,Estrgbs,rcongbs] = fn_fcstidcnd(valuecon,stepcon,varcon,nstepsm,...
                     nconstr,eq_ms,nvar,lags,phil,0,1,yforehat,imf3shat,A0hat,Bhat,...
                     nfqm,TLindx,TLnumber,nCms,eq_Cms);
            ycforegbse = [fdates ycforegbs];
         else
            %=== Unconditional forecasts.
            if nexo<2
               yforegbs = fn_forecastsim(Bhat,A0hat,phil,[nvar lags nfqm]);
                              % unconditional with simulated shocks; nfqm-by-nvar, in log
            else
               yforegbs = fn_forecastsim(Bhat,A0hat,phil,[nvar lags nfqm],nexo,Xfexoe(:,3:end));
                              % unconditional with simulated shocks; nfqm-by-nvar, in log
            end
            ycforegbse = [fdates yforegbs];
         end

         yafgbse = [yact1e; ycforegbse];  % actual and unconditional forecast.  Need actual to compute annual growth rates
         %===== Converted to mg (or qg) and yg
         [yafyrggbse,jnk,yafqmyggbse] = fn_datana(yafgbse,q_m,vlistlog(1:nlogeno),vlistper(1:npereno));
                           % actual and unconditional forecast growth rates

         yforeyrggbse = fn_dataext([fdates(1,1) 0],[fdates(end,1) 0],yafyrggbse);
         yforeqmyggbse = fn_dataext(fdates(1,:),fdates(end,:),yafqmyggbse);
                                   % unconditional forecast growth rates with the dates
         if Foreband
            temp = yforeyrggbse(:,3:end);
            yrgu1 = yrgu1+temp(:);     % 1st moment  annual growth rate in unconditional forecast
            yrgu2 = yrgu2+temp(:).^2;  % 2nd moment  annual growth rate in unconditional forecast
            temp = yforeqmyggbse(:,3:end);
            qmyrgu1 = qmyrgu1+temp(:);  % 1st moment annualized (to q or m in the last year) growth rate in unconditional forecast
            qmyrgu2 = qmyrgu2+temp(:).^2;  % 2nd moment annualized (to q or m in the last year) growth rate in unconditional forecast
         end
      end
      %
   end

   if IndxFore
      if Foreband
         yrgu1 = yrgu1/ndraws1;
         yrgu2 = yrgu2/ndraws1;
         yrgustd = abs(sqrt(yrgu2 - yrgu1.^2));   % standard deviation.  Abs is used to gauard agaist almost zero but imaginary numbers.
         yrgustd(find(yrgustd<sqrt(eps)))=1e-04;   % set to a very small number for zero standard deviations.
         yrguRange=zeros(nfyr*nvar,2);  % 2: min and max
         yrguRange(:,1) = yrguhat(:)-nstd*abs(yrgustd);
         yrguRange(:,2) = yrguhat(:)+nstd*abs(yrgustd);
         yrguhbin = (yrguRange(:,2) - yrguRange(:,1)) ./ ninv;
         %
         qmyrgu1 = qmyrgu1/ndraws1;
         qmyrgu2 = qmyrgu2/ndraws1;
         qmyrgustd = abs(sqrt(qmyrgu2 - qmyrgu1.^2));   % standard deviation.  Abs is used to gauard agaist almost zero but imaginary numbers.
         qmyrgustd(find(qmyrgustd<sqrt(eps)))=1e-04;   % set to a very small number for zero standard deviations.
         qmyrguRange=zeros(nfqm*nvar,2);  % 2: min and max
         qmyrguRange(:,1) = qmyrguhat(:)-nstd*abs(qmyrgustd);
         qmyrguRange(:,2) = qmyrguhat(:)+nstd*abs(qmyrgustd);
         qmyrguhbin = (qmyrguRange(:,2) - qmyrguRange(:,1)) ./ ninv;
      end
   end


   %*** Second set of draws of forecasts to get the distributions
   tic
   for draws = 1:ndraws2
      if ~mod(draws,nbuffer)
         draws
         %  fwriteid = fopen('outA0.bin','a');
         %  count = fwrite(fwriteid,A0hatw,'double');
         %  status = fclose('all');
      end

      if IndxFore   % draws of out-of-sample forecasts
         if nconstr %WZ-DSL conditional forecasts.
            [ycforegbs,Estrgbs,rcongbs] = fn_fcstidcnd(valuecon,stepcon,varcon,nstepsm,...
                     nconstr,eq_ms,nvar,lags,phil,0,1,yforehat,imf3shat,A0hat,Bhat,...
                     nfqm,TLindx,TLnumber,nCms,eq_Cms);
            ycforegbse = [fdates ycforegbs];
         else
            %=== Unconditional forecasts.
            if nexo<2
               yforegbs = fn_forecastsim(Bhat,A0hat,phil,[nvar lags nfqm]);
                              % unconditional with simulated shocks; nfqm-by-nvar, in log
            else
               yforegbs = fn_forecastsim(Bhat,A0hat,phil,[nvar lags nfqm],nexo,Xfexoe(:,3:end));
                              % unconditional with simulated shocks; nfqm-by-nvar, in log
            end
            ycforegbse = [fdates yforegbs];
         end



         yafgbse = [yact1e; ycforegbse];  % actual and conditional forecast.  Need actual to compute annual growth rates
         %===== Converted to mg (or qg) and yg
         [yafyrggbse,jnk,yafqmyggbse,jnk2,yafyqgbse, jnk3, yaflevelgbse] = fn_datana2(yafgbse,q_m,vlistlog(1:nlogeno),vlistper(1:npereno));
                              % actual and conditional forecast growth rates
         if (draws==1)
            yaflevelgbse_draws = zeros(size(yaflevelgbse,1),size(yaflevelgbse,2),ndraws2);
            yafyqgbse_draws = zeros(size(yafyqgbse,1),size(yafyqgbse,2),ndraws2);
         end
         yaflevelgbse_draws(:,:,draws) = yaflevelgbse;
         yafyqgbse_draws(:,:,draws) = yafyqgbse;


         yforeyrggbse = fn_dataext([fdates(1,1) 0],[fdates(end,1) 0],yafyrggbse);
         yforeqmyggbse = fn_dataext(fdates(1,:),fdates(end,:),yafqmyggbse);
                                   % unconditional forecast growth rates with the dates
         if Foreband
            yrgugbs = yforeyrggbse(:,3:end);
            yrgumean = yrgumean+yrgugbs;     % mean.   nfyr-by-nvar
            yrgucnt = fn_empdfsort(yrgucnt,yrgugbs,yrguRange(:,1),yrguhbin,ninv);
                          % count or distribution.   ninv+2-by-nfyr*nvar
            %
            qmyrgugbs = yforeqmyggbse(:,3:end);
            qmyrgumean = qmyrgumean+qmyrgugbs;     % mean.  nfqm-by-nvar
            qmyrgucnt = fn_empdfsort(qmyrgucnt,qmyrgugbs,qmyrguRange(:,1),qmyrguhbin,ninv);
                          % count or distribution.  ninv+2-by-nfqm*nvar
         end
      end
      %
   end
   timeminutes=toc/60

   if IndxFore
      if Foreband
         yrgumean=yrgumean/ndraws2;
         qmyrgumean=qmyrgumean/ndraws2;
      end
   end
   %
end

save outdata_distcnt n0 np nvar imstp nfyr nfqm Indxcol ndraws1 ndraws2 timeminutes forelabel...
                 q_m fdates yafyrghate yafqmyghate yact2yrge yact2qmyge ...
                 nCsk ninv IndxParR IndxOvR IndxAp IndxImf IndxFore ...
                 Aband bhat bmean bstd bcnt bhbin bRange ...
                 Apband ghat gmean gstd gcnt ghbin gRange ...
                 Imfband imfhat imfmean imfstd imfcnt imfhbin imfRange ...
                 Foreband yrguhat yrgumean yrgustd yrgucnt yrguhbin yrguRange ...
                 qmyrguhat qmyrgumean qmyrgustd qmyrgucnt qmyrguhbin qmyrguRange ...
                 imfgbs_draws

%======= The following is for debugging =============
%yrgu1=reshape(yrgu1,nfyr,nvar);
%yrgustd = reshape(yrgustd,nfyr,nvar);
%[yrgu1(:,4) yrgumean(:,4) yrguhat(:,4) yrgustd(:,4)]
%%
%qmyrgu1=reshape(qmyrgu1,nfqm,nvar);
%qmyrgustd = reshape(qmyrgustd,nfqm,nvar);
%[qmyrgu1(:,4) qmyrgumean(:,4) qmyrguhat(:,4) qmyrgustd(:,4)]

%======= Error bands of forecasts.
if IndxFore   % draws of out-of-sample forecasts
   %--- Level
   yaflevelgbse_draws_sort = sort(yaflevelgbse_draws,3);
   yaflevelgbse_low70 = yaflevelgbse_draws_sort(:,:,0.15*ndraws2);
   yaflevelgbse_high70 = yaflevelgbse_draws_sort(:,:,0.85*ndraws2);
   save outdataout_fore_level_low70.prn yaflevelgbse_low70 -ascii
   save outdataout_fore_level_high70.prn yaflevelgbse_high70 -ascii
   %--- Annualized quarterly growth rate.
   yafyqgbse_draws_sort = sort(yafyqgbse_draws,3);
   yafyqgbse_low70 = yafyqgbse_draws_sort(:,:,0.15*ndraws2);
   yafyqgbse_high70 = yafyqgbse_draws_sort(:,:,0.85*ndraws2);
   save outdataout_fore_qg_low70.prn yafyqgbse_low70 -ascii
   save outdataout_fore_qg_high70.prn yafyqgbse_high70 -ascii
end



%------------------------
%-- Error bands or densities of impulse responses
%------------------------
%if (IndxImf & Imfband)
%   ncom = 10;   % The number of combined intervals.  A large ncom implies a large size of the bin.
%   %            % Must set ncom so that ninv can be divided
%   kIndx=[1 100 imstp*nvar^2];  % index for selected variables.  Length(kIndx)<=size(bcnt,2).
%   lval = [];   % length(kIndx)-element vector of the lowest values on the axis;
%         % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
%   hval = [];   % length(kIndx)-element vector of the highest values on the axis;
%         % The vector corresponds to kIndx.  This option is disenabled by setting it to [].
%   [imfpdf,imfpo,imfprob] = fn_histpdfcnt(imfcnt,ndraws2,imfRange(:,1),imfhbin,ninv,0,ncom,kIndx,lval,hval);
%            % ninv+2-by-imstp*nvar^2 matrix of the densities (histograms) of variables
%   [imferrl,imferrh,imferrl1,imferrh1] = fn_demarw(imfpo,imfprob);
%            % .90 and .68 error bands
%   %*** plot impulse responses
%   imferrl = reshape(imferrl,imstp,nvar^2);
%   imferrh = reshape(imferrh,imstp,nvar^2);
%   imferrl1 = reshape(imferrl1,imstp,nvar^2);
%   imferrh1 = reshape(imferrh1,imstp,nvar^2);
%   figure
%   fn_imcerrgraph(imfhat,imferrl1,imferrh1,nvar,imstp,xlab,ylab,1,[6 12 24 36])  % .68 bands.
%   subtitle('.68 Error Bands')
%   figure
%   fn_imcerrgraph(imfhat,imferrl,imferrh,nvar,imstp,xlab,ylab,1,[6 12 24 36])  % .90 bands.
%   subtitle('.90 Error Bands')
%   %

%   %====== <<>> Plotting only responses to spending and tax shocks. ======
%   if (indx_model==0)
%      indx_sel_shks = [1 2];
%      indx_sel_vars = [1:5];
%      indx_rescale_vars = [1 2];
%   elseif (indx_model==1)
%      indx_sel_shks = [1];
%      indx_sel_vars = [1:2];
%      indx_rescale_vars = [1];
%   end
%   irfs_expand_all = zeros(imstp, nvar, nvar, 3); %nvar variables, nvar shocks, 3: hat, low band, and high band.
%   irfs_expand_all(:,:,:,1) = imf3hat;
%   irfs_expand_all(:,:,:,2) = reshape(imferrl1,size(imferrl1,1),nvar,nvar);
%   irfs_expand_all(:,:,:,3) = reshape(imferrh1,size(imferrh1,1),nvar,nvar);
%   if indx_rescale
%      irfs_expand_all(:,:,indx_sel_shks(1),1) = irfs_expand_all(:,:,indx_sel_shks(1),1)*0.01 / irfs_expand_all(1,indx_rescale_vars(1),indx_sel_shks(1),1);
%      irfs_expand_all(:,:,indx_sel_shks(2),1) = -irfs_expand_all(:,:,indx_sel_shks(2),1)*0.01 / irfs_expand_all(1,indx_rescale_vars(2),indx_sel_shks(2),1);
%      %
%      irfs_expand_all(:,:,indx_sel_shks(1),2) = irfs_expand_all(:,:,indx_sel_shks(1),2)*0.01 / irfs_expand_all(1,indx_rescale_vars(1),indx_sel_shks(1),2);
%      irfs_expand_all(:,:,indx_sel_shks(2),2) = -irfs_expand_all(:,:,indx_sel_shks(2),2)*0.01 / irfs_expand_all(1,indx_rescale_vars(2),indx_sel_shks(2),2);
%      %
%      irfs_expand_all(:,:,indx_sel_shks(1),3) = irfs_expand_all(:,:,indx_sel_shks(1),3)*0.01 / irfs_expand_all(1,indx_rescale_vars(1),indx_sel_shks(1),3);
%      irfs_expand_all(:,:,indx_sel_shks(2),3) = -irfs_expand_all(:,:,indx_sel_shks(2),3)*0.01 / irfs_expand_all(1,indx_rescale_vars(2),indx_sel_shks(2),3);
%   end
%   figure('PaperPosition',[0.25 1 8 9]); %[left bottom witdh height]  Tip: bottom + height = 10.0
%   scaleout = fn_multigraphn_ver2(irfs_expand_all(:,indx_sel_vars, indx_sel_shks,:),xlab(indx_sel_shks),ylab(indx_sel_vars),'Quarter',[4 8 16,24],0,0,length(indx_sel_vars), length(indx_sel_shks));
%   print -depsc2 outfig_irfs.eps
%end
%

if (IndxImf)
   imf_demarcated = prctile(imfgbs_draws,[5 16 50 84 95],1);
   imferrl = imf_demarcated(1,:)';
   imferrl1 = imf_demarcated(2,:)';
   imfmedian = imf_demarcated(3,:)';
   imferrh1 = imf_demarcated(4,:)';
   imferrh = imf_demarcated(5,:)';
   %*** plot impulse responses
   imferrl = reshape(imferrl,imstp,nvar^2);
   imferrh = reshape(imferrh,imstp,nvar^2);
   imfmedian = reshape(imfmedian,imstp,nvar^2);
   imferrl1 = reshape(imferrl1,imstp,nvar^2);
   imferrh1 = reshape(imferrh1,imstp,nvar^2);
   if indx_model==2 || indx_model==4;
       figure
       fn_imcerrgraph_lastshk(imfmedian,imferrl1,imferrh1,nvar,imstp,xlab,ylab,1,[q_m 2*q_m 3*q_m 4*q_m],1)  % .68 bands.
       subtitle('.68 Error Bands')   
       figure
       fn_imcerrgraph_lastshk(imfmedian,imferrl,imferrh,nvar,imstp,xlab,ylab,1,[q_m 2*q_m 3*q_m 4*q_m],1)  % .90 bands.
       subtitle('.90 Error Bands')

       %<<>> Printing out ad hoc figure
       figure('PaperPosition',[0.25 1 8 9]); %[left bottom witdh height]  Tip: bottom + height = 10.0
       fn_imcerrgraph_lastshk(imfmedian,imferrl1,imferrh1,nvar,imstp,xlab,ylab,1,[q_m 2*q_m 3*q_m 4*q_m],1)  % .68 bands.
       %subtitle('.68 Error Bands')   
       print -depsc2 outfig_irfs_qu_68bands.eps 
       print -dpdf outfig_irfs_qu_68bands.pdf 
       close   
   else
       figure
       fn_imcerrgraph_1shk(imfmedian,imferrl1,imferrh1,nvar,imstp,xlab,ylab,1,[q_m 2*q_m 3*q_m 4*q_m],1)  % .68 bands.
       subtitle('.68 Error Bands')   
       figure
       fn_imcerrgraph_1shk(imfmedian,imferrl,imferrh,nvar,imstp,xlab,ylab,1,[q_m 2*q_m 3*q_m 4*q_m],1)  % .90 bands.
       subtitle('.90 Error Bands')
       if imfmedian(1,1)<0;
           imfmed = -imfmedian(:,1:4);
           imf_errl = -imferrh(:,1:4);
           imf_errh = -imferrl(:,1:4);
       else
           imfmed = imfmedian(:,1:4);
           imf_errl = imferrl(:,1:4);
           imf_errh = imferrh(:,1:4);
       end;
       xlswrite('irfs_out',[imfmed, imf_errl, imf_errh]);

       %<<>> Printing out ad hoc figure
       figure('PaperPosition',[0.25 1 8 9]); %[left bottom witdh height]  Tip: bottom + height = 10.0
       fn_imcerrgraph_1shk(imfmedian,imferrl1,imferrh1,nvar,imstp,xlab,ylab,1,[q_m 2*q_m 3*q_m 4*q_m],1)  % .68 bands.
       %subtitle('.68 Error Bands')   
       print -depsc2 outfig_irfs_qu_68bands.eps 
       print -dpdf outfig_irfs_qu_68bands.pdf 
       close 
       
   end;
   

%%====== <<>> Plotting only responses to spending and tax shocks. ======
%if (indx_model==0)
%   indx_sel_shks = [1 2];
%   indx_sel_vars = [1:5];
%   indx_rescale_vars = [1 2];
%elseif (indx_model==1)
%   indx_sel_shks = [1];
%   indx_sel_vars = [1:2];
%   indx_rescale_vars = [1];
%elseif (indx_model==2)
%   indx_sel_shks = [1 2];
%   indx_sel_vars = [1:3];
%   indx_rescale_vars = [];   
%lseif (indx_model==3)
%  indx_sel_shks = [1 2];
%  indx_sel_vars = [1:3];
%  indx_rescale_vars = [];
%else
%   error('Must specify the model type indx_model');
%end
%irfs_expand_all = zeros(imstp, nvar, nvar, 3); %nvar variables, nvar shocks, 3: hat, low band, and high band.
%if (0)
%   irfs_expand_all(:,:,:,1) = imf3hat;
%else
%   irfs_expand_all(:,:,:,1) = reshape(imfmedian,size(imfmedian,1),nvar,nvar);
%end
%irfs_expand_all(:,:,:,2) = reshape(imferrl1,size(imferrl1,1),nvar,nvar);
%irfs_expand_all(:,:,:,3) = reshape(imferrh1,size(imferrh1,1),nvar,nvar);
%if indx_rescale
%   irfs_expand_all(:,:,indx_sel_shks(1),1) = irfs_expand_all(:,:,indx_sel_shks(1),1)*0.01 / irfs_expand_all(1,indx_rescale_vars(1),indx_sel_shks(1),1);
%   irfs_expand_all(:,:,indx_sel_shks(2),1) = -irfs_expand_all(:,:,indx_sel_shks(2),1)*0.01 / irfs_expand_all(1,indx_rescale_vars(2),indx_sel_shks(2),1);
%   %
%   irfs_expand_all(:,:,indx_sel_shks(1),2) = irfs_expand_all(:,:,indx_sel_shks(1),2)*0.01 / irfs_expand_all(1,indx_rescale_vars(1),indx_sel_shks(1),2);
%   irfs_expand_all(:,:,indx_sel_shks(2),2) = -irfs_expand_all(:,:,indx_sel_shks(2),2)*0.01 / irfs_expand_all(1,indx_rescale_vars(2),indx_sel_shks(2),2);
%   %
%   irfs_expand_all(:,:,indx_sel_shks(1),3) = irfs_expand_all(:,:,indx_sel_shks(1),3)*0.01 / irfs_expand_all(1,indx_rescale_vars(1),indx_sel_shks(1),3);
%   irfs_expand_all(:,:,indx_sel_shks(2),3) = -irfs_expand_all(:,:,indx_sel_shks(2),3)*0.01 / irfs_expand_all(1,indx_rescale_vars(2),indx_sel_shks(2),3);
%end
%figure('PaperPosition',[0.25 1 8 9]); %[left bottom witdh height]  Tip: bottom + height = 10.0
%scaleout = fn_multigraphn_ver2(irfs_expand_all(:,indx_sel_vars, indx_sel_shks,:),xlab(indx_sel_shks),ylab(indx_sel_vars),'Quarter',[4 8 16,24],0,0,length(indx_sel_vars), length(indx_sel_shks));
%subtitle('Median and .68 bands')
%print -depsc2 outfig_irfs.eps
end


