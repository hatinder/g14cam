F2B2RK2T1000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM2B2RK2THETA1000.txt';
F2BR2K2T2000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM2B2RK2THETA2000.txt';

delimiterIn = ' ';
headerlinesIn = 1;
RK2T21000 = importdata(F2B2RK2T1000, delimiterIn,headerlinesIn);
RK2T22000 = importdata(F2BR2K2T2000, delimiterIn,headerlinesIn);
x = [ sin(RK2T21000.data(:,2)),  sin(RK2T21000.data(:,2))+sin(RK2T22000.data(:,2))];
y = [-cos(RK2T21000.data(:,2)), -cos(RK2T21000.data(:,2))-cos(RK2T22000.data(:,2))];
% Convert radians to degrees
ang = [RK2T21000.data(:,2),RK2T22000.data(:,2)]*180/pi;
N=2000;

figure;
subplot(2,1,1);
xlabel('time (sec)'); ylabel('angle (\circ)');
T=RK2T21000.data(:,1);
for i=1:N
   
    subplot(2,1,1);
    plot(T,ang, 'LineWidth', 2);
    line(T(i), ang(i,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
    line(T(i), ang(i,2), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
    xlabel('time (sec)'); ylabel('angle (deg)');

   % The bottom plot shows the animation of the double pendulum
   subplot(2,1,2);
   plot([0, x(i,1);x(i,1), x(i,2)], [0, y(i,1);y(i,1), y(i,2)], '.-', 'MarkerSize', 20, 'LineWidth', 2);
   axis equal; axis([-2 2 -2 2]);
   title(sprintf('Time: %0.2f sec', T(i)));
   %disp(i);
   drawnow;
end
 

% hold on;
% 
% plot(RK2T1000.data(:,2));
% plot(log10(FIMLP000.data(:,1)),log10(abs(FIMLP000.data(:,2))));
% grid on;
% hold off;
%legend('\phi_0(x)','location','best');
% xlabel('n');
% ylabel('log_{10}( | I-I_n | )');
% title('log_{10}( | I-I_n | ) plot against n'); 
