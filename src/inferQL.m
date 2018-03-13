% inferQL.m : This script calculates the sufficient statistics for the
% variational posterior over factor loading matrix entries.
%
% modified by
% Emanuele Sansone GCNU 15/12/14

n = size(Y,2);
p = size(Y,1);
s = size(Lm,2);

for t = 1:s
  kt = size(Lm{t},2);
  mean_Lambda = [mean_mcl zeros(p,kt-1)];
  num{t} = [nu_mcl repmat(a./b{t},[p 1]) ];
  temp = Xm{t}.*repmat(sum(Qns(:,t,:),3)',kt,1);
  T2 = Xcov{t}*sum(sum(Qns(:,t,:),3),1) + Xm{t}*temp';
  T3 = diag(psii)* Y *temp';
  for q = 1:p
    Lcov{t}(:,:,q) = inv( diag(num{t}(q,:)) + psii(q)*T2);
    Lm{t}(q,:) = ( T3(q,:) + mean_Lambda(q,:)*diag(num{t}(q,:)) ) * Lcov{t}(:,:,q);
  end
end


