% circle.m : The function generates points uniformly 
% from a circle centered at (x1,y1) with radius rc
%
% Added by
% Emanuele Sansone GCNU 15/12/14


function x=circle(x1,y1,rc)
    a=2*pi*rand;
    r=sqrt(rand);
    x(1)=(rc*r)*cos(a)+x1;
    x(2)=(rc*r)*sin(a)+y1;
end