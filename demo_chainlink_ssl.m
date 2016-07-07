% Added by
% Emanuele Sansone GCNU 15/12/14

close all; clear all; clc;
disp('Loading data (chainlink distribution)');
load data/Chainlink.mat;
chainlink = chainlink';
[chainlink ppp] = preprocess(chainlink);
disp('Running Variational Bayesian Algorithm');
figure;
plot_ssl_data(chainlink,chainlink_labels);
title('Ground trouth')
[chainlink chainlink_labels] = sample_dataset(chainlink,chainlink_labels,200);
net = vbmfa(chainlink,chainlink_labels,3,0,1,10);
plot_output_data(chainlink,net);
title('Predicted labels')