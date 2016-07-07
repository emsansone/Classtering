% inferQbeta.m : This script infers the sufficient statistics for the
% variational posterior for the Dirichlet parameter $\beta$.
%
% Added by
% Emanuele Sansone GCNU 15/12/14

s = size(Lm,2);
pg = gamma/K*ones(s,K);

g = pg + reshape(sum(Qns,1),[s K]);