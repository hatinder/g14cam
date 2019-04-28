m = 1;         % mass
L = 1;         % link length
theta1 = pi/6; % initial angle for theta 1 (in rad)
theta2 = pi/6; % initial angle for theta 2 (in rad)
t = linspace(0, 10, 300);  % simulate for 10 seconds with 300 points

% Solving ODE of a double pendulum
[T,Y] = ode45(@(t, x) double_pendulum(t, x, m, L), t, [theta1, theta2, 0, 0]);

% Calculating joint coordinates for animation purposes
x = [ L*sin(Y(:,1)),  L*sin(Y(:,1))+L*sin(Y(:,2))];
y = [-L*cos(Y(:,1)), -L*cos(Y(:,1))-L*cos(Y(:,2))];

% Convert radians to degrees
ang = Y(:,1:2)*180/pi;

figure;
subplot(2,1,1);
xlabel('time (sec)'); ylabel('angle (\circ)');
tic;    % start timing
for id = 1:length(T)
   % The top plot shows a time series of link angles
   subplot(2,1,1);
   plot(T,ang, 'LineWidth', 2);
   line(T(id), ang(id,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
   line(T(id), ang(id,2), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
   xlabel('time (sec)'); ylabel('angle (deg)');

   % The bottom plot shows the animation of the double pendulum
   subplot(2,1,2);
   plot([0, x(id,1);x(id,1), x(id,2)], [0, y(id,1);y(id,1), y(id,2)], ...
      '.-', 'MarkerSize', 20, 'LineWidth', 2);
   axis equal; axis([-2*L 2*L -2*L 2*L]);
   title(sprintf('Time: %0.2f sec', T(id)));

   drawnow;
end
fprintf('Animation (Regular): %0.2f sec\n', toc);