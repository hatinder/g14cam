PRB3AMATRIX000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\MATRIX064.txt';
PRB3AVECTOR000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\VECTOR064.txt';

delimiterIn = ' ';
% headerlinesIn = 1;
DPRB3AMATRIX000 = importdata(PRB3AMATRIX000, delimiterIn);
DPRB3AVECTOR000 = importdata(PRB3AVECTOR000, delimiterIn);

x = DPRB3AVECTOR000; 
y = DPRB3AVECTOR000;
% disp(x);
[X,Y] = meshgrid(x,y);
Z=DPRB3AMATRIX000;
contour(X,Y,Z,32);
mesh(X,Y,Z);
% disp(meshgrid(x,y));
 