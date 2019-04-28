FRK2LP000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM2ARK2ERROR000.txt';
FIMLP000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM2AIMERROR000.txt';

delimiterIn = ' ';
headerlinesIn = 1;
RK2LP000 = importdata(FRK2LP000, delimiterIn,headerlinesIn);
FIMLP000 = importdata(FIMLP000, delimiterIn,headerlinesIn);

hold on;
plot(log10(RK2LP000.data(:,1)),log10(abs(RK2LP000.data(:,2))));
plot(log10(FIMLP000.data(:,1)),log10(abs(FIMLP000.data(:,2))));
grid on;
hold off;
%legend('\phi_0(x)','location','best');
% xlabel('n');
% ylabel('log_{10}( | I-I_n | )');
% title('log_{10}( | I-I_n | ) plot against n'); 
