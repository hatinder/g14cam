LP000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM2AIMERROR000.txt';

delimiterIn = ' ';
headerlinesIn = 1;
LP000 = importdata(LP000, delimiterIn,headerlinesIn);

plot(log10(LP000.data(:,1)),log10(LP000.data(:,2)));
grid on;
%legend('\phi_0(x)','location','best');
xlabel('n');
ylabel('log_{10}( | I-I_n | )');
title('log_{10}( | I-I_n | ) plot against n'); 