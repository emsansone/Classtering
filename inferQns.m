% inferQns.m : This script calculates the variational posterior for
% the class-conditional responsibilities for the data. It also removes
% from the model any components which have zero (or below threshold)
% total responsibility.
%
% modified by
% Emanuele Sansone GCNU 15/12/14

n = size(Y,2);
p = size(Y,1);
s = size(Lm,2); % number of components
Qns_old = Qns;

allocprobs = sum(sum(Qns,2),3); % vector obtained by summing over S and K
col = 0; logQns = [];
for t = 1:s
  col = col + 1;
  kt = size(Lm{t},2);
  LmpsiiLm = Lm{t}'*diag(psii)*Lm{t};
  temp = LmpsiiLm + reshape(reshape(Lcov{t},kt*kt,p)*psii,kt,kt);

  logQns(:,col) = -.5*( +sum(Y.*(diag(psii)*(Y-2*Lm{t}*Xm{t})),1)' ... % This term belongs to trace (psii)
      +reshape(temp,kt*kt,1)'*reshape(Xcov{t},kt*kt,1) ...             % This term belongs to trace (psii)
      +sum( Xm{t}.*(temp*Xm{t}) ,1)' ...                               % This term belongs to trace (psii)
      +trace(Xcov{t}(2:end,2:end)) ...
      +sum( Xm{t}(2:end,:).*Xm{t}(2:end,:) ,1)' ...
      -2*sum(log(diag(chol(Xcov{t}(2:end,2:end))))) ...                % Computing the determinant term (Chapter 2 MacKay)
      );
end %t
  
logQns = logQns + ones(n,1)*digamma(u);

% Adding new terms
logQns3 = repmat(logQns,1,1,K);
tmp = reshape(g,[1 s*K]);
tmp = digamma(tmp);
tmp2 = sum(g,2);
tmp2 = repmat(digamma(tmp2'),[1 K]);
g_temp = reshape(tmp,[1,s,K])-reshape(tmp2,[1,s,K]);
class = [log(reshape(double(Y_labels),N,1,K)); zeros(n-N,1,K)];
logQns3 = logQns3 + repmat(g_temp,n,1,1) + repmat(class,1,s,1);
logQns3 = logQns3 - repmat(max(max(logQns,[],2),[],3),1,s,K); % Force all numbers to be negative
Qns = exp(logQns3);
Qns = Qns .* (  repmat((allocprobs./sum(sum(Qns,2),3)),1,s,K)  ); % scaling with new values of Qns


% check for any empty components, and remove them if allowed
if removal==1
  [dummy,empty_t] = find(sum(sum(Qns,1),3) < 1 );
  num_died = size(empty_t,2);
  % if there exist components to remove, do so
  if num_died > 0
    remain = 1:size(Lm,2); remain(empty_t) = []; 
    fprintf('\nRemoving component: '); fprintf('%2i, ',empty_t);
    % ascertain if a cophd
    num_children = size(intersect(empty_t,[parent size(Lm,2)]),2);
    if num_children~=num_died
      cophd=0; % i.e. require reordering
    else
      if num_children == 1
        fprintf('\nChild of parent has died')
        cophd = 1;
        decleft = candorder(1:pos-1)>parent;
        decright = candorder(pos+1:end)>parent;
        if empty_t == parent
          fprintf('\nDead component is original parent - need to change ordering');
          fprintf('\nOld order'); fprintf(' %i',candorder);
          decrement = candorder(1:pos-1)>parent;
          candorder = [candorder(1:pos-1)-decleft size(Lm,2)-1 candorder(pos+1:end)-decright];
          fprintf('\nNew order'); fprintf(' %i',candorder);
        else % empty_t == size(Lm,2)
          fprintf('\nDead component is newly born - no change to component ordering');
        end
      end
      if num_children == 2
        cophd = 1;
        decleft = candorder(1:pos-1)>parent;
        decright = candorder(pos+1:end)>parent;
        fprintf('\nBoth children died - need to change ordering');
        fprintf('\nOld order'); fprintf(' %i',candorder);
        candorder = [candorder(1:pos-1)-decleft candorder(pos+1:end)-decright];
        fprintf('\nNew order'); fprintf(' %i',candorder);
        pos = pos-1; % because 2 components died.
      end
    end
    
    Lm = Lm(1,remain);
    Lcov = Lcov(1,remain);
    Xm = Xm(1,remain);
    Xcov = Xcov(1,remain);
    b = b(1,remain);
    u = u(remain);
    pu = pu(remain);
    g = g(remain,:);
    pg = pg(remain,:);
    s = size(remain,2);
    Qns = Qns(:,remain,:);
    %inferQns % N.B. this may cause stacking.
  end
end

if size(Qns_old,2) == size(Qns,2)
  dQns = abs(Qns-Qns_old);
  dQns_sagit = (sum(dQns,1)./sum(Qns,1)); % The percentage absolute movement
else
  dQns_sagit = ones(1,size(Qns,2),size(Qns,3)); % i.e. 100% agitated

end

Qns = Qns./repmat(sum(sum(Qns,2),3),[1 size(Qns,2) size(Qns,3)]); % Normalisation
Ps = sum(sum(Qns,1),3)/size(Qns,1);
Pk = sum(sum(Qns,1),2)/size(Qns,1);
