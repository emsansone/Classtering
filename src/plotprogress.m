% Added by
% Emanuele Sansone GCNU 15/12/14

%fig = figure(3);

%subplot('position',[.7 .05 .28 .9]);
%bar(Ps);
%set(1,'doublebuffer','on');
%title('Mixing Proportions')

%orbh = subplot('position',[0.02 .05 .6 .9]);
orbh = subplot('position',[0.02 .05 .98 .9]);

switch p
  case 2
    color = {'*r','*g','*c','*m','*y','*k',...
             'or','og','oc','om','oy','ok'};
    for i = 1:K
        idx = Y_labels(:,i)'>0;
        plot(Y(1,idx),Y(2,idx),color{i},'MarkerSize',10);
        hold on;
    end
    if N < n
        plot(Y(1,(N+1):end),Y(2,(N+1):end),'.k','MarkerSize',10);
    end
  case 3
    color = {'*r','*g','*c','*m','*y','*k',...
             'or','og','oc','om','oy','ok'};
    for i = 1:K
        idx = Y_labels(:,i)'>0;
        plot3(Y(1,idx),Y(2,idx),Y(3,idx),color{i},'MarkerSize',10);
        hold on;
    end
    if N < n
        plot3(Y(1,(N+1):end),Y(2,(N+1):end),Y(3,(N+1):end),'.k','MarkerSize',10);
    end
end
hold on;
h1 = zeros(1,size(Lm,2));
for t = 1:size(Lm,2);
  h1(t) = plot_gaussian(Lm{t}(:,2:end)*Lm{t}(:,2:end)'+diag(1./psii),Lm{t}(:,1),t,15); hold on 
end
%title('Classtering')
axis equal
%axis off
set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
set(gca, 'ZTickLabelMode', 'manual', 'ZTickLabel', []);
grid on;
hold off;

%view(orbh,375,38);
drawnow