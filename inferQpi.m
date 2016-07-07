% inferQpi.m : This script infers the sufficient statistics for the
% variational posterior for the Dirichlet parameter $\bpi$.
%
% modified by
% Emanuele Sansone GCNU 15/12/14

s = size(Lm,2);
pu = alpha/s *ones(1,s);

u = pu + sum(sum(Qns,1),3);


