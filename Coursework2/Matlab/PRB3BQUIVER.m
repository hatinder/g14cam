PRB3BMATU000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PRB3BMATRIXU000.txt';
PRB3BMATV000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PRB3BMATRIXV000.txt';

delimiterIn = ' ';
% headerlinesIn = 1;
DPRB3BMATU000 = importdata(PRB3BMATU000, delimiterIn);
DPRB3BMATV000 = importdata(PRB3BMATV000, delimiterIn);

u = DPRB3BMATU000; 
v = DPRB3BMATV000;
N=33;
x = linspace(0,1,N);
y = linspace(0,1,N);
% disp(x);

[X,Y] = meshgrid(x,y);
hold on;
h2=quiver(X,Y,u,v);
set(h2,'AutoScale','on', 'AutoScaleFactor',2)
axis equal square
% colormap hsv;
hold off;
% Z=DPRB3AMATRIX000;

% contour(X,Y,Z,32);
% mesh(X,Y,v);
% disp(meshgrid(x,y));