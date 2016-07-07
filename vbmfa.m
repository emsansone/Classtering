%Classtering: Variational Mixture of Factor Analysers model 
%
%net = vbmfa(Y,Y_labels,maxdim,pcaflag,Fflag,dsp,net);
%
% Y - p x n matrix of (normalised) observations (see preprocess.m)
% Y_labels - N x K matrix of labels (where N <= n)
% maxdim - maximum factor dimensionality (default p-1)
% pcaflag - noise model - FA:0 (default), PCA:1
% Fflag - calculate and display F at each iteration - off:0 (default), on:1
% dsp - display interval (default never)
% net - if specified hparams will be extracted, and, if present,
%  params will be extracted as well.
%
%net is a structure, consisting of structures:
% 
% hparams - hyperparameters over the model
%            [mean_mcl nu_mcl],psii,pa,pb,alpha,gamma
% params  - inferred parameters of the mixture components
%            Lm,Lcov,u,g
% hidden  - inferred hidden states
%            Xm,Xcov,Qns,a,b
%
% Fhist   - history of F, very verbose if Fflag set to 1
%
% Emanuele Sansone GCNU 15/12/14
% Initial version (Matthew J. Beal GCNU 15/02/01)

function net = vbmfa(Y,Y_labels,maxdim,pcaflag,Fflag,dsp,net);

[N K] = size(Y_labels);
[p n] = size(Y);

strng = './Atom/Atom';

if nargin<3
  k = (p-1) + 1; % i.e. that sufficient to capture covariance, plus one for the bias.
else
  k = maxdim + 1;
end
if nargin<4, pcaflag = 0; end
if nargin<5, Fflag = 0; end
if nargin<6, dsp = pi; end
if nargin==7 % given hparams, also perhaps params
  Fhist = net.Fhist;
  mean_mcl = net.hparams.mcl(:,1);
  nu_mcl = net.hparams.mcl(:,2);
  psii = net.hparams.psii;
  pa = net.hparams.pa;
  pb = net.hparams.pb;
  alpha = net.hparams.alpha;
  gamma = net.hparams.gamma;
  if ~isempty(net.params) % given params
    Lm = net.params.Lm;
    Lcov = net.params.Lcov;
    u = net.params.u;
    g = net.params.g;
  else % generate random params
    k = (p-1) +1;
    Lm{1} = randn(p,k); Lcov{1} = repmat(eye(k),[1 1 p]); u = 1;
    g = 1/K*ones(1,K);
  end
  net.hidden = [];
else % nothing given, generate random everything
  net = struct( ...
      'model','variational MFA', ...
      'hparams',struct('mcl',[],'psii',[],'pa',[],'pb',[],'alpha',[],'gamma',[]), ...
      'params',struct('Lm',[],'Lcov',[],'u',[],'g',[]), ...
      'hidden',struct('Xm',[],'Xcov',[],'Qns',[],'a',[],'b',[]), ...
      'Fhist',[]);
  pa = 1; pb = 1; alpha = 1; gamma = 1; psii = ones(p,1);
  for i = 1:K
      Lm{i}(:,1) = mean(Y(:,find(Y_labels(:,i)==1)),2);
      matrix = chol(posdef_matrix(cov(Y(:,find(Y_labels(:,i)==1))')-diag(psii)))';
      Lm{i}(:,2:k) = matrix(:,1:k-1);
      Lcov{i} = repmat(eye(k),[1 1 p]);
  end
  mean_mcl = mean(Y,2); nu_mcl = (1./std(Y,0,2)).^2;
  u = 1/K*ones(1,K); g = 1/K*ones(K,K);

end

s = size(Lm,2); % number of components
pu = alpha*ones(1,s)/s; % uniform probability over each component
pg = gamma*ones(s,K)/K;

if (p==2 || p==3) && dsp ~= 0 && K <= 6
    figure(2);
    plot_ssl_data(Y,Y_labels);
    title('Training dataset');
    Ps = pu;
    fprintf('\nINITIALIZATION\n');
end


% Tweakables and initialisations
balance = 1;
psimin = 1e-5;
maintol = 6; % the greater these numbers,
finetol = 8; %  the finer learning
birthtol = 1e-2;
maxtries = 1; tries = 0;
removal = 1; % if on, insists removal of components that have close to zero responsibility
F = -inf;
Fhist = [];
it = 1; I = eye(k-1,k-1);

% Infer hidden variables and calculate initial F
tol = exp(-maintol);

% Initializing Qns
Qns = ones(n,s,K)/(s*K);
for i = 1:K
    matrix = (ones(1,s,K) - ones(1,s,K)*balance)/((s*K-1)*(1-balance)+balance);
    Qns(find(Y_labels(:,i)==1),:,:) = repmat(matrix,[size(find(Y_labels(:,i)==1),1),1,1]);
    Qns(find(Y_labels(:,i)==1),i,i) = balance; 
end

inferQX, inferQns, inferQnu, learn


% Learning based on the inizialized configuration
it = 1; epoch = 0;
Fcalc; Ftarg = F; dF = 0;
dQns_sagit = tol;
while ~isempty(find(dQns_sagit>=exp(-maintol)))
    learn
    if Fflag == 1
      dF_old = dF; Fcalc;
      fprintf('\nit:%4i F:%4.2f dF:%4.2f ',it,F,dF);
    end
    display = 1;
    if dsp > it
        display = rem(dsp,it);
    else
        display = rem(it,dsp);
    end
    if (p==2 || p==3) && display == 0 && K <= 6
        plotprogress;
    end
end

workspace = {psii Lm Lcov Qns a b u g};

candorder = [1:size(Lm,2) 0]; pos = 1; parent = candorder(pos);


while parent ~= 0;
  
  epoch = epoch+1; cophd = 0;
  Fcalc; Ftarg = F; dF = 0;
  F_birth = F;
  dobirth

  tol = maintol;
  dQns_sagit = tol;
  it = 1;

  % Learning for new configuration
  while ~isempty(find(dQns_sagit>=exp(-tol)))
    learn
    if Fflag == 1
      dF_old = dF; Fcalc;
      fprintf('\nit:%4i F:%4.2f dF:%4.2f ',it,F,dF);
    end
    display = 1;
    if dsp > it
        display = rem(dsp,it);
    else
        display = rem(it,dsp);
    end
    if (p==2 || p==3) && display == 0 && K <= 6
        plotprogress;
    end
  end
  
  % Analyse whether to accept Epoch
  Fcalc;
    
  if F > Ftarg
    fprintf('\nACCEPTING F= %4.2f (%4.2f)',F,Ftarg);
    % If child of parent has died do not reorder components
    if cophd == 1
      fprintf('\nChild of parent has died, no reordering.');
      pos = pos+1;
    else
      ordercands % sets "candorder", with last entry as a ZERO (termination)
      tries = 0; % resets the tries back to zero as well as the pos.
    end
    workspace = {psii Lm Lcov Qns a b u g};
    inferQX; infermcl; learn
    Fcalc; Ftarg = F; if exist('orbh')&(p==3); orbit(orbh,360,30); end
  else
    fprintf('\nREJECTING F= %4.2f (%4.2f) ',F,Ftarg);
    [psii,Lm,Lcov,Qns,a,b,u,g] = deal(workspace{:});
    inferQX; infermcl; learn
    pos = pos +1;
  end
  
  fprintf('\nit:%4i Completion: %2.1f%%',it,( tries/maxtries+(pos-1)/(size(candorder,2)-1)/maxtries )*100);
  fprintf('\n');
  
  % Move along the candidate list; if we reach the end start over.
  parent = candorder(pos);
  if parent == 0;
    tries = tries+1;
    if tries ~= maxtries
      fprintf('\nEnd of ordering reached, reordering and trying more splits\n')
      inferQX; inferQns; 
      Fcalc; ordercands;
      parent = candorder(pos);
    end
  end

end

fprintf('\nOptimisation complete with %2i components\n\n',size(Lm,2));
if (p==2 || p==3) && dsp ~= 0 && K <= 6
    plotprogress; 
    string = sprintf([strng, '%d_%d.pdf'],epoch,it);
end



net.hparams.mcl = [mean_mcl nu_mcl];
net.hparams.psii = psii;
net.hparams.pa = pa;
net.hparams.pb = pb;
net.hparams.alpha = alpha;
net.hparams.gamma = gamma;
net.params.Lm = Lm;
net.params.Lcov = Lcov;
net.params.u = u;
net.params.g = g;
net.hidden.Xm = Xm;
net.hidden.Xcov = Xcov;
net.hidden.Qns = Qns;
net.hidden.a = a;
net.hidden.b = b;
net.Fhist = Fhist;