% Added by
% Emanuele Sansone GCNU 15/12/14

function [sigma] = posdef_matrix(sigma) 

epsilon = 1e-6;
zero = 1e-8;

[~, err] = cholcov(sigma, 0);

if (err ~= 0) 
    [v d] = eig(sigma);
    d=diag(d);
    d( d <= zero ) = epsilon;
    d=diag(d);
    sigma = v*d*v';
end