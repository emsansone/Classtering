% dobirth.m : This script performs a birth operation. It has
% considerable flexibility in the way it does this, along with a few
% tweakables, but here the operation is to split the
% *responsibilities* of the parent component.

% Please email m.beal@gatsby.ucl.ac.uk for more information.
% modified by
% Emanuele Sansone GCNU 15/12/14

s = size(Lm,2);
n = size(Y,2);
k = size(Lm{parent},2);

s = s+1;
t = s; % t is the identity of the newborn.

hardness = 1; % 1 means hard, 0.5 means soft, 0 is equiv to 1 (opposite direction)

fprintf('\nCreating comp:%2i, as hybrid of comp:%2i.',t,parent);

Lm{t} = Lm{parent};
[U,S,D] = svd(Lm{parent}(:,2:end)*Lm{parent}(:,2:end)'+diag(1./psii));
delta_vector = D(:,1)*S(1,1);

% Selection of samples which contribute to the parent component
Qns_temp = squeeze(sum(Qns,3));
Qns_temp(find(Qns_temp >= 0.5)) = 1;
Qns_temp(find(Qns_temp < 0.5)) = 0;
index = find(Qns_temp(:,parent)==1);

% Qns birth
assign = sign( D(:,1)'*(Y-repmat(Lm{parent}(:,1),1,n)) ); % size 1 x n
pos_ind = find(assign == 1);
neg_ind = find(assign == -1);
% Reassign those one side of vector to the child, t,
% whilst the rest remain untouched. Positive are sent to child
Qns(pos_ind',t,:) = hardness*Qns(pos_ind',parent,:);
Qns(pos_ind',parent,:) = (1-hardness)*Qns(pos_ind',parent,:);
Qns(neg_ind',t,:) = (1-hardness)*Qns(neg_ind',parent,:);
Qns(neg_ind',parent,:) = hardness*Qns(neg_ind',parent,:);

% Set all features of t to those of t_parent
assign = sign( D(:,1)'*(Y(:,index)-repmat(Lm{parent}(:,1),1,length(index))) ); % size 1 x n
pos_ind = find(assign == 1);
neg_ind = find(assign == -1);
Y_temp = Y(:,index);
Y_parent = (Y_temp(:,neg_ind') - repmat(Lm{parent}(:,1),1,length(neg_ind)))./repmat(sqrt(psii),1,length(neg_ind));
[COEFF, SCORE] = princomp(Y_parent');
Lm{parent}(:,2:end) = COEFF(:,1:k-1);
Y_t = (Y_temp(:,pos_ind') - repmat(Lm{parent}(:,1),1,length(pos_ind)))./repmat(sqrt(psii),1,length(pos_ind));
[COEFF, SCORE] = princomp(Y_t');
Lm{t}(:,2:end) = COEFF(:,1:k-1);

Lcov{t} = Lcov{parent};
Lm{t}(:,1)      = Lm{parent}(:,1) + delta_vector;
Lm{parent}(:,1) = Lm{parent}(:,1) - delta_vector;
b{t} = b{parent};
u(parent) = u(parent)/2; u(t) = u(parent);
pu = alpha/s *ones(1,s);
g = [g;g(parent,:)];
pg = gamma/K*ones(t,K);

% Update the model
infermcl, inferQX
Fcalc; 
dF = abs(F-F_birth);
while dF > birthtol;
    learn_without_Qns
    F_birth = F; Fcalc; dF = abs(F-F_birth);
end


allcomps
cophd = 0;

