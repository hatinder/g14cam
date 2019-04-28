PRB2BIM1TH000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PRB2BIM1TH000.txt';
PRB2BIM2TH000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PRB2BIM2TH000.txt';

delimiterIn = ' ';
headerlinesIn = 1;
DPRB2BIM1TH000 = importdata(PRB2BIM1TH000, delimiterIn,headerlinesIn);
DPRB2BIM2TH000 = importdata(PRB2BIM2TH000, delimiterIn,headerlinesIn);

 hold on;
% 
plot(DPRB2BIM1TH000.data(:,1),DPRB2BIM1TH000.data(:,2));
plot(DPRB2BIM2TH000.data(:,1),DPRB2BIM2TH000.data(:,2));
%plot(log10(FIMLP000.data(:,1)),log10(abs(FIMLP000.data(:,2))));
grid on;
hold off;
legend('First Initial Condition','Second Initial Condition','location','best');
xlabel('Time');
ylabel('Hamiltonian');
title('Implicit Midpoint'); 
