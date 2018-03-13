% demo_gauss.m: This is the demo for data distributed according to 5
% Gaussians in a semi-supervised setting
%
% Added by
% Emanuele Sansone GCNU 15/12/14

close all; clear all; clc

addpath(genpath(pwd))

disp('Generating 200 samples from 5 Gaussians');
mu = [-2 20; 6 13; 10 5; -10 0; -20 -2];
%mu = [1 2; 4 1; 1 5; -1 0; -2 -2];
sigma(:,:,1) = [1 0.4; 0.4 1];
sigma(:,:,2) = [0.9 -0.8; -0.8 0.9];
sigma(:,:,3) = [2 0.5; 0.5 0.9];
sigma(:,:,4) = [1 0.6; 0.6 0.9];
sigma(:,:,5) = [1 0.9; 0.9 1.5];
%rng default;
obj = gmdistribution(mu,sigma);
Y = random(obj, 600);
Y = Y';
idx = Y(1,:)>-2;
Y_labels(:,1) = idx';
Y_labels(:,2) = (~idx)';
disp('Running Variational Bayesian Algorithm');
[Y,ppp] = preprocess(Y);
%plot(Y(1,:),Y(2,:),'o');
figure;
plot_ssl_data(Y,Y_labels);
title('Ground trouth');
[Y Y_labels] = sample_dataset(Y,Y_labels,200);
net = vbmfa(Y,Y_labels,1,0,5,5);
plot_output_data(Y,net);
title('Predicted labels');