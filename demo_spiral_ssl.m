% Added by
% Emanuele Sansone GCNU 15/12/14

close all; clear all; clc;
disp('Loading data (spiral distribution)');
load data/Spiral.mat;
spiral = spiral';
[spiral ppp] = preprocess(spiral);
disp('Running Variational Bayesian Algorithm');
figure;
plot_ssl_data(spiral,spiral_labels);
title('Ground truth');
[spiral spiral_labels] = sample_dataset(spiral,spiral_labels,700);
net = vbmfa(spiral,spiral_labels,3,0,1,10);
plot_output_data(spiral,net);
title('Predicted labels');