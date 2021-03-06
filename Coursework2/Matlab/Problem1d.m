LP000 = '..\Code\cmake-build-debug\PROBLEM1D000.txt';

delimiterIn = ' ';
headerlinesIn = 1;
LP000 = importdata(LP000, delimiterIn,headerlinesIn);

plot(LP000.data(:,1),log10(LP000.data(:,2)));
grid on;
%legend('\phi_0(x)','location','best');
xlabel('n');
ylabel('log_{10}( | I-I_n | )');
title('log_{10}( | I-I_n | ) plot against n'); 