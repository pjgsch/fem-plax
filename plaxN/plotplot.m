%**********************************************************************

function plotplot(xx,yy,r1,r2,xl,yl,tt);

nplots = size(xx,1);
%line(1).style = '-';
%line(2).style = ':';
%line(3).style = '--';
%line(4).style = '-.';
line(1).color = 'b';
line(2).color = 'k';
line(3).color = 'r';
line(4).color = 'm';

if     nargin==2, r1=1; r2=size(xx,2); xl=''; yl=''; tt=''; 
elseif nargin==3, r2=size(xx,2); xl=''; yl=''; tt=''; 
elseif nargin==4, xl=''; yl=''; tt=''; 
elseif nargin==5, yl=''; tt=''; 
elseif nargin==6, tt=''; 
end;

if r2<0, r2 = size(xx,2); end;

for np=1:nplots
%  plot( xx(np,r1:r2),yy(np,r1:r2) );
  plot( xx(np,r1:r2),yy(np,r1:r2),line(np).color );
%  plot( xx(np,r1:r2),yy(np,r1:r2),line(np).style );
  hold on;
end;
hold off;
grid on;

xlabel(xl); ylabel(yl); title(tt);

%**********************************************************************
