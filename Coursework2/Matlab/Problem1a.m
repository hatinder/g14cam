LP000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM1A000.txt';
LP001 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM1A001.txt';
LP002 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM1A002.txt';
LP003 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM1A003.txt';
LP004 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM1A004.txt';

delimiterIn = ' ';
headerlinesIn = 1;
LP000 = importdata(LP000, delimiterIn,headerlinesIn);
LP001 = importdata(LP001, delimiterIn,headerlinesIn);
LP002 = importdata(LP002, delimiterIn,headerlinesIn);
LP003 = importdata(LP003, delimiterIn,headerlinesIn);
LP004 = importdata(LP004, delimiterIn,headerlinesIn);

plot(LP000.data(:,1),LP000.data(:,2),LP001.data(:,1),LP001.data(:,2),LP002.data(:,1),LP002.data(:,2),LP003.data(:,1),LP003.data(:,2),LP004.data(:,1),LP004.data(:,2))

legend('\phi_0(x)','\phi_1(x)','\phi_2(x)','\phi_3(x)','\phi_4(x)','location','bestoutside');
xlabel('x')
ylabel('\phi_n(x)')
title('Legendre Polynomial')