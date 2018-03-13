% This function takes as input a test data set, hyperparameters and
% parameters of a model, and outputs the inferred hidden state
% posterior as well as a calculation of the lower bound for the
% evidence for the test data.
%
% ensure data is preprocessed using same 'ppparams' as was used on
% training set.
%
% [tehidden,F] = infer(Ytest,net);
%
% modified by
% Emanuele Sansone GCNU 15/12/14

function [tehidden,Ftest] = infer(Y,net);

psii = net.hparams.psii;
pa = net.hparams.pa;
pb = net.hparams.pb;
alpha = net.hparams.alpha;
gamma = net.hparams.gamma;
Lm = net.params.Lm;
Lcov = net.params.Lcov;
u = net.params.u;
g = net.params.g;

s = size(Lm,2);
n = size(Y,2); p = size(Y,1);
K = size(g,2);
allcomps
pa_mcl =pa; pb_mcl = 25*p;  pu = alpha*ones(1,s)/s; pg = gamma*ones(1,K)/K;
tol = 6; dQns_sagit_tol = exp(-tol); removal = 0;
Y_labels = ones(n,K);
N = n;

% calculate P for the test data. We can't find marginal for P
% analytically.  Content to report few terms of \cF as a lower bound
% for the log predictive density.

Qns = ones(n,s,K)/(s*K);
inferQnu
allcomps
% No need to iteratively update posteriors QX and Qns, as QX posterior
% is not a function of Qns posterior.
inferQX
inferQns

% calculate last few terms of $\mathcal{F}$, i.e. terms 2, 3 and 4
PsK = reshape(sum(Qns,1),s,K)/n;
Fmatrix3 = zeros(4,s,K,n);
tempdig_u = digamma(u);
tmp = reshape(g,[1 s*K]);
tempdig_g = reshape(digamma(tmp),[s,K]);
tempdig_usum = digamma(sum(u,2));
tempdig_usum_g = digamma(sum(g,2)')';
Qnsmod = Qns; Qnsmod(find(Qnsmod==0)) = 1;
Y_labelsmod = Y_labels; Y_labelsmod(find(Y_labelsmod==0)) = 1;


for t = 1:s
  kt = size(Lm{t},2);
  for j = 1:K
       
       Fmatrix3(1,t,j,:) = Qns(:,t,j)'.*( -log(Qnsmod(:,t,j)') ...
                                         +ones(n,1,1)'*(tempdig_u(t)-tempdig_usum) ...
                                         +ones(n,1,1)'*(tempdig_g(t,j)-tempdig_usum_g(t)));
  
       Fmatrix3(2,t,j,:) = +.5*(kt-1)*PsK(t,j) ...
                           -.5*trace(Xcov{t}(2:end,2:end))*Qns(:,t,j)' -.5*sum(Xm{t}(2:end,:).^2,1).*Qns(:,t,j)' ...
                           +.5*PsK(t,j)*( log(det(Xcov{t}(2:end,2:end))) );     

       temppp = Lm{t}'*diag(psii)*Lm{t} + reshape(reshape(Lcov{t},kt*kt,p)*psii,kt,kt);             
       Fmatrix3(3,t,j,:) = -.5*PsK(t,j)*( -log(det(diag(psii))) +p*log(2*pi) ) ...
              -.5*trace( Xcov{t} *temppp )*Qns(:,t,j)' ...
              -.5*Qns(:,t,j)'.*sum(Xm{t}.*(temppp'*Xm{t}),1) ...
              -.5*sum((diag(psii)*(Y.*repmat(Qns(:,t,j)',p,1))).*(Y-2*Lm{t}*Xm{t}),1);
       
       Fmatrix3(4,t,j,:) = Qns(:,t,j)'.*log(double(Y_labelsmod(:,j)')); % $\log(c_n|I_n)$
  end
end

Ftest = squeeze(sum(sum(sum(Fmatrix3,1),2),3))';

tehidden.Xm = Xm;
tehidden.Xcov = Xcov;
tehidden.Qns = Qns;
tehidden.a = a;
tehidden.b = b;