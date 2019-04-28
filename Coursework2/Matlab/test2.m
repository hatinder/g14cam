figure;
subplot(2,1,1);
plot(T, ang, 'LineWidth', 2);
hh1(1) = line(T(1), ang(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(2) = line(T(1), ang(1,2), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
xlabel('time (sec)'); ylabel('angle (deg)');

subplot(2,1,2);
hh2 = plot([0, x(1,1);x(1,1), x(1,2)], [0, y(1,1);y(1,1), y(1,2)], ...
      '.-', 'MarkerSize', 20, 'LineWidth', 2);
axis equal
axis([-2*L 2*L -2*L 2*L]);
ht = title(sprintf('Time: %0.2f sec', T(1)));

tic;     % start timing
for id = 1:length(T)
   % Update XData and YData
   set(hh1(1), 'XData', T(id)          , 'YData', ang(id, 1));
   set(hh1(2), 'XData', T(id)          , 'YData', ang(id, 2));
   set(hh2(1), 'XData', [0, x(id, 1)]  , 'YData', [0, y(id, 1)]);
   set(hh2(2), 'XData', x(id, :)       , 'YData', y(id, :));
   set(ht, 'String', sprintf('Time: %0.2f sec', T(id)));

   drawnow;
end
fprintf('Animation (Smart update): %0.2f sec\n', toc);