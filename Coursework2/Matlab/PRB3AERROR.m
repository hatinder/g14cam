PRB3AERR000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PRB3AERROR000.txt';

delimiterIn = ' ';
headerlinesIn = 1;
DPRB3AERR000 = importdata(PRB3AERR000, delimiterIn,headerlinesIn);

 hold on;
% 
plot(log(DPRB3AERR000.data(:,1)),log(DPRB3AERR000.data(:,2))); 
%plot(log10(FIMLP000.data(:,1)),log10(abs(FIMLP000.data(:,2))));
grid on;
hold off;
%legend('First Initial Condition','Second Initial Condition','location','best');
ylabel('log(|u-\hat u|)');
xlabel('N');
title('Problem 3A'); 
