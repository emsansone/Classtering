% plot_output_data.m: This function plots data with the estimated labels
%
%   Y - p x n matrix containing n samples described by p-dimensional
%       feature vectors
%   net - structure obtained as a result from the algorithm
%
% Added by
% Emanuele Sansone GCNU 15/12/14
%

function plot_output_data(Y, net)

figure;
labels = sum(net.hidden.Qns,2);
labels = squeeze(labels);
labels(find(labels<.5))=0;
labels(find(labels>=.5))=1;
labels = logical(labels);
plot_ssl_data(Y,labels);