F2BRK2T1000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM2BRK2THETA1000.txt';
F2BRK2T2000 = 'E:\edubitbucket\g14cam\Coursework2\Code\cmake-build-debug\PROBLEM2BRK2THETA2000.txt';
delimiterIn = ' ';
headerlinesIn = 1;
RK2T1000 = importdata(F2BRK2T1000, delimiterIn,headerlinesIn);
RK2T2000 = importdata(F2BRK2T2000, delimiterIn,headerlinesIn);
T=RK2T1000.data(:,1);
%N=size(T);
N=500;
x = [ sin(RK2T1000.data(:,2)),  sin(RK2T1000.data(:,2))+sin(RK2T2000.data(:,2))];
y = [-cos(RK2T1000.data(:,2)), -cos(RK2T1000.data(:,2))-cos(RK2T2000.data(:,2))];
hh2 = plot([0, x(1,1);x(1,1), x(1,2)], [0, y(1,1);y(1,1), y(1,2)], '.-', 'MarkerSize', 20, 'LineWidth', 2);
axis equal; axis([-2*L 2*L -2*L 2*L]);
ht = title(sprintf('Time: %0.2f sec', T(1)));

writerObj = VideoWriter('2BRK2_i.mp4','MPEG-4'); % Name it.
writerObj.FrameRate = 60; % How many frames per second.
writerObj.Quality=80;
open(writerObj);
for i=1:N
   set(hh2(1), 'XData', [0, x(i, 1)]  , 'YData', [0, y(i, 1)]);
   set(hh2(2), 'XData', x(i, :)       , 'YData', y(i, :));
   set(ht, 'String', sprintf('Time: %0.2f sec', T(i)));
   frame=getframe(gcf);
   writeVideo(writerObj,frame);
end
close(writerObj);

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
