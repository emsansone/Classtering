% plot_ssl_data.m: This function plots data
%
%   Y - p x n matrix containing n samples described by p-dimensional
%       feature vectors
%   labels - N x K matrix of labels (K classes)
%
% N can be lower than n
%
% Added by
% Emanuele Sansone GCNU 15/12/14
%

function plot_ssl_data(Y,labels)

[N K] = size(labels);
[p n] = size(Y);

orbh = subplot('position',[0.02 .05 .98 .9]);

switch p
  case 2
    color = {'*r','*g','*c','*m','*y'};
    for i = 1:K
        idx = labels(:,i)'>0;
        plot(Y(1,idx),Y(2,idx),color{i},'MarkerSize',10);
        hold on;
    end
    if N < n
        plot(Y(1,(N+1):end),Y(2,(N+1):end),'.k','MarkerSize',10);
    end
  case 3
    color = {'*r','*g','*c','*m','*y'};
    for i = 1:K
        idx = labels(:,i)'>0;
        plot3(Y(1,idx),Y(2,idx),Y(3,idx),color{i},'MarkerSize',10);
        hold on;
    end
    if N < n
        plot3(Y(1,(N+1):end),Y(2,(N+1):end),Y(3,(N+1):end),'.k','MarkerSize',10);
    end    
end
%title('Classtering')
axis equal
%axis off
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
set(gca, 'ZTickLabelMode', 'manual', 'ZTickLabel', []);
grid on;
hold off;
