% Added by
% Emanuele Sansone GCNU 15/12/14

close all; clear all; clc;

addpath(genpath(pwd))

disp('Loading data (ring distribution)');
load data/Ring.mat;
ring = ring';
[ring ppp] = preprocess(ring);
disp('Running Variational Bayesian Algorithm');
figure;
plot_ssl_data(ring,ring_labels);
title('Ground truth');
[ring ring_labels] = sample_dataset(ring,ring_labels,40);
net = vbmfa(ring,ring_labels,2,0,1,10);
plot_output_data(ring,net);
title('Predicted labels');