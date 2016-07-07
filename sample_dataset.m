% sample_dataset.m: This function generates a dataset of only M*K labeled
% samples
%
%   Y - p x n matrix containing n samples described by p-dimensional
%       feature vectors
%   labels - n x K matrix of labels (K classes)
%   M - number of samples with labels for each class
%
% N can be lower than n
%
% Added by
% Emanuele Sansone GCNU 15/12/14
%

function [out_Y, out_labels] = sample_dataset(Y, labels, M)

[n K] = size(labels);

out_labels = [];
out_Y = [];
indices = [];
for i = 1:K
    idx = find(labels(:,i)==1)';
    idy = [];
    for j = 1:M
        tmp = randi(length(idx));
        idy = [idy, idx(tmp)];
        idy = unique(sort(idy));
        if length(idy) < j
            idx(tmp) = [];
            idy = [idy, idx(randi(length(idx)))];
            idy = unique(sort(idy));            
        end
    end
    y = Y(:,idy);
    out_Y = [out_Y, y];
    out_labels = [out_labels; labels(idy,:)];
    indices = [indices, idy];
end
indices = sort(indices);
Y(:,indices) = [];
out_Y = [out_Y, Y];