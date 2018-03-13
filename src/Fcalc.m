% Fcalc.m : This script calculates the lower bound on the log evidence
% using the formula for $\cF$.
%
% modified by
% Emanuele Sansone GCNU 15/12/14

s = size(Lm,2);
p = size(Y,1);

PsK = reshape(sum(Qns,1),s,K)/n;
F_old = F;

Fmatrix = zeros(2,s);
Fmatrix3 = zeros(4,s,K);
tempdig_a = digamma(a);
tempdig_u = digamma(u);
tmp = reshape(g,[1 s*K]);
tempdig_g = reshape(digamma(tmp),[s,K]);
tempdig_usum = digamma(sum(u,2));
tempdig_usum_g = digamma(sum(g,2)')';
pu = alpha/s *ones(1,s);
Qnsmod = Qns; Qnsmod(find(Qnsmod==0)) = 1;
Y_labelsmod = Y_labels; Y_labelsmod(find(Y_labelsmod==0)) = 1;

FmatrixKLpi = -kldirichlet(u,pu);
FmatrixKLbeta = zeros(s,1);


for t = 1:s
  kt = size(Lm{t},2);
  for j = 1:K
       temp_alt = Xm{t}.*repmat(Qns(:,t,j)',kt,1);
       Xcor_t_Qns_weighted = Xcov{t}*sum(Qns(:,t,j),1)  +  Xm{t}*temp_alt';
  
       Fmatrix3(1,t,j) = +sum( Qns(:,t,j).*( -log(Qnsmod(:,t,j)) ...
                                             +ones(n,1,1)*(tempdig_u(t)-tempdig_usum) ...
                                             +ones(n,1,1)*(tempdig_g(t,j)-tempdig_usum_g(t))) );
  
       Fmatrix3(2,t,j) = +.5*n*(kt-1)*PsK(t,j) ...
                      -.5*trace( Xcor_t_Qns_weighted(2:end,2:end) ) ...
                      +.5*n*PsK(t,j)*( ...
                        2*sum(log(diag(chol(Xcov{t}(2:end,2:end))))) );
      
       Fmatrix3(3,t,j) = -.5*n*PsK(t,j)*( -log(det(diag(psii))) +p*log(2*pi) ) ...
                      -.5*trace( Lm{t}'*diag(psii)*Lm{t} * Xcor_t_Qns_weighted ) ...
                      -.5*trace( reshape(reshape(Lcov{t},kt*kt,p)*psii,kt,kt) * Xcor_t_Qns_weighted ) ...
                      -.5*trace( diag(psii)* ((ones(p,1)*Qns(:,t,j)').*Y) * (Y-2*Lm{t}*Xm{t})' );
       class = [log(double(Y_labelsmod(:,j))); zeros(n-N,1)];
       Fmatrix3(4,t,j) = sum(Qns(:,t,j).*class); % $\log(c_n|I_n)$
  end
  Fmatrix(2,t) = -klgamma([a],[b{t}],pa,pb);
  FmatrixKLbeta(t) = -kldirichlet(g(t,:),pg(t,:));
end

% Fmatrix(1,:) reimplementation
for t = 1:s
  f1 = 0;
  priorlnnum{t} = [log(nu_mcl) repmat(tempdig_a-log(b{t}),[p 1])];
  priornum{t} = [nu_mcl repmat(a./b{t},[p 1])];
  f1 = f1 + sum(sum(priorlnnum{t}));
  for q = 1:p
    f1 = f1 + log(det(Lcov{t}(:,:,q))) - size(Lcov{t}(:,:,q),1);
    f1 = f1 - ( diag(Lcov{t}(:,:,q))'+Lm{t}(q,:).^2 ) * priornum{t}(q,:)';
    f1 = f1 - priornum{t}(q,1)*( -2*Lm{t}(q,1)*mean_mcl(q,1) + mean_mcl(q,1).^2);
  end
  Fmatrix(1,t) = f1/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = sum(sum(sum(Fmatrix3))) + sum(sum(Fmatrix)) + FmatrixKLpi + sum(FmatrixKLbeta);
dF = F-F_old; 
Fhist = [Fhist; it F];
