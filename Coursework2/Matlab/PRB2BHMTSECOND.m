PRB2BHMTRK2000 = '..\Code\cmake-build-debug\PRB2BRK21TH000.txt';
PRB2BHMTIM000 = '..\Code\cmake-build-debug\PRB2BIM1TH000.txt';

delimiterIn = ' ';
headerlinesIn = 1;
DPRB2BHMTRK2000 = importdata(PRB2BHMTRK2000, delimiterIn,headerlinesIn);
DPRB2BHMTIM000 = importdata(PRB2BHMTIM000, delimiterIn,headerlinesIn);

 hold on;
% 
plot(DPRB2BHMTRK2000.data(:,1),DPRB2BHMTRK2000.data(:,2));
plot(DPRB2BHMTRK2000.data(:,1),DPRB2BHMTIM000.data(:,2));
%plot(log10(FIMLP000.data(:,1)),log10(abs(FIMLP000.data(:,2))));
grid on;
hold off;
legend('Runge Kutta 2','Implicit Midpoint','location','best');
xlabel('Time');
ylabel('Hamiltonian');
title('First Initial Condition'); 
