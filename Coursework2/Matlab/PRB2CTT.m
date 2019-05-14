PRB2BIM1TH000 = '..\Code\cmake-build-debug\PRB2CTT000.txt';

delimiterIn = ' ';
headerlinesIn = 1;
DPRB2BIM1TH000 = importdata(PRB2BIM1TH000, delimiterIn,headerlinesIn);

 hold on;
% 
plot(DPRB2BIM1TH000.data(:,1)*180/pi,DPRB2BIM1TH000.data(:,2));
%plot(log10(FIMLP000.data(:,1)),log10(abs(FIMLP000.data(:,2))));
grid on;
hold off;
%legend('First Initial Condition','Second Initial Condition','location','best');
ylabel('Time');
xlabel('\Theta_2');
title('Problem 2C'); 
