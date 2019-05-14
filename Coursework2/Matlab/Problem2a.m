FAU000 = '..\Code\cmake-build-debug\PROBLEM2ARK2U000.txt';
FAU001 = '..\Code\cmake-build-debug\PROBLEM2ARK2U001.txt';
FAU002 = '..\Code\cmake-build-debug\PROBLEM2ARK2U002.txt';
FAU003 = '..\Code\cmake-build-debug\PROBLEM2ARK2U003.txt';
FAU004 = '..\Code\cmake-build-debug\PROBLEM2ARK2U004.txt';
FAU005 = '..\Code\cmake-build-debug\PROBLEM2ARK2U005.txt';

delimiterIn = ' ';
headerlinesIn = 1;
AU000 = importdata(FAU000, delimiterIn,headerlinesIn);
AU001 = importdata(FAU001, delimiterIn,headerlinesIn);
AU002 = importdata(FAU002, delimiterIn,headerlinesIn);
AU003 = importdata(FAU003, delimiterIn,headerlinesIn);
AU004 = importdata(FAU004, delimiterIn,headerlinesIn);
AU005 = importdata(FAU005, delimiterIn,headerlinesIn);

hold on;
plot(AU000.data(:,1),AU000.data(:,2));
plot(AU001.data(:,1),AU001.data(:,2));
plot(AU002.data(:,1),AU002.data(:,2));
plot(AU003.data(:,1),AU003.data(:,2));
plot(AU004.data(:,1),AU004.data(:,2));
plot(AU005.data(:,1),AU005.data(:,2));
grid on;
hold off;
%legend('\phi_0(x)','location','best');
% xlabel('n');
% ylabel('log_{10}( | I-I_n | )');
% title('log_{10}( | I-I_n | ) plot against n'); 