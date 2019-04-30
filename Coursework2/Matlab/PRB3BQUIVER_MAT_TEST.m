PRB3BMATC000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PRB3BMATRIXC000.txt';
PRB3BMATF000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PRB3BMATRIXF000.txt';

delimiterIn = ' ';
% headerlinesIn = 1;
DPRB3BMATC000 = importdata(PRB3BMATC000, delimiterIn);
DPRB3BMATF000 = importdata(PRB3BMATF000, delimiterIn);
% DPRB3BMATV000 = importdata(PRB3BMATV000, delimiterIn);
C=sparse(DPRB3BMATC000);
F=DPRB3BMATF000;

% disp(C);
[l,u,p]=lu(C);
u1=l\p*F;
ufinal=u\u1;
disp(ufinal);
% u = DPRB3BMATU000; 
% v = DPRB3BMATV000;
% N=7;
% x = linspace(0,1,N);
% y = linspace(0,1,N);
% % disp(x);
% 
% [X,Y] = meshgrid(x,y);
% hold on;
% h2=quiver(X,Y,u,v);
% set(h2,'AutoScale','on', 'AutoScaleFactor',2)
% axis equal square
% % colormap hsv;
% hold off;
% % Z=DPRB3AMATRIX000;
% 
% % contour(X,Y,Z,32);
% % mesh(X,Y,v);
% disp(meshgrid(x,y));