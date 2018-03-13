% demo_cake_ssl.m: This is the demo for data distributed according to a
% cake in a semi-supervised scenario
%
% Added by
% Emanuele Sansone GCNU 15/12/14

close all; clear all; clc

addpath(genpath(pwd))

x = [];
MEANX = 2;
MEANY = 1;
%rng default;
for i = 1:1000
x = [x circle(MEANX,MEANY,3)'];
end
x = x';
if 0
RADIUS = 2;
x1 = @(t) t.*cos(t);
y1 = @(t) t.*sin(t);

x2 = @(t) RADIUS.*t.*cos(t);
y2 = @(t) RADIUS.*t.*sin(t);

x3 = @(t) -t.*cos(t);
y3 = @(t) -t.*sin(t);

x4 = @(t) -RADIUS.*t.*cos(t);
y4 = @(t) -RADIUS.*t.*sin(t);
end

id1 = x(:,1)<MEANX & x(:,2)>MEANY;
id2 = x(:,1)>MEANX & x(:,2)<MEANY;
Y_labels(:,1) = (id1 | id2);
Y_labels(:,2) = (~(id1 | id2));
Y = x';
[Y,ppp] = preprocess(Y);
clear x id1 id2 i MEANX MEANY ppp;
figure;
plot_ssl_data(Y,Y_labels);
title('Ground trouth')
[Y Y_labels] = sample_dataset(Y,Y_labels,50);
net = vbmfa(Y,Y_labels,2,0,1,10);
plot_output_data(Y,net);
title('Predicted labels')