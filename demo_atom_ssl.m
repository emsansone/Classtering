% Added by
% Emanuele Sansone GCNU 15/12/14

close all; clear all; clc;

addpath(genpath(pwd))

disp('Loading data (atom distribution)');
load data/Atom.mat;
atom = atom';
[atom ppp] = preprocess(atom);
disp('Running Variational Bayesian Algorithm');
figure;
plot_ssl_data(atom,atom_labels);
title('Ground trouth')
[atom atom_labels] = sample_dataset(atom,atom_labels,40);
net = vbmfa(atom,atom_labels,3,0,1,10);
plot_output_data(atom,net);
title('Predicted labels')