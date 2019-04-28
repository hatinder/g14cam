% x = linspace(0,1);
% y = linspace(0,1);
% [X,Y] = meshgrid(x,y);
% Z =cos(4*pi*X)*cos(4*pi*Y);
% contour(X,Y,Z)
N=6;
x = linspace(0,1,N);
y = linspace(0,1,N);
% [X,Y]=meshgrid(x,y);
% disp(X);
% disp(Y);
% surf(x,y,cos(4*pi*x).*cos(4*pi*y));
% disp(x);
[X,Y] = meshgrid(x,y);
% disp(X);disp(Y);
Z=zeros(N,N);
for i=1:N
  for j=1:N
    Z(j,i) = cos(4*pi*x(j))*cos(4*pi*x(i));
%     temp=sprintf("(%d,%d) - X:%f, Y:%f , Z:%f",j,i,x(j),x(i),Z(j,i));
%     disp(temp);
  end
end
mesh(X,Y,Z);
% disp(Z);
contour(X,Y,Z,64);