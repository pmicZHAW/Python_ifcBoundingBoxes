clc, clear variables

% R = C_GB = Rz * Ry * Rx

rpy = [-30, 45, 60]*pi/180;
R = calcRn([0;0;1], rpy(3)) * ...
    calcRn([0;1;0], rpy(2)) * ...
    calcRn([1;0;0], rpy(1));

p = [1, 2, 3].'

fprintf('%0.8f %0.8f %0.8f\n', R.')

corners = [[-1 -0.5 -0.25]; ...
           [-1 -0.5  0.25]; ...
           [-1  0.5 -0.25]; ...
           [-1  0.5  0.25]; ...
           [ 1 -0.5 -0.25]; ...
           [ 1 -0.5  0.25]; ...
           [ 1  0.5 -0.25]; ...
           [ 1  0.5  0.25]]
       
unique(corners(:,[2 3]), 'rows')
unique(corners(:,[1 3]), 'rows')
unique(corners(:,[1 2]), 'rows')

t_corners = corners * R.' + p.'

unique(t_corners(:,[2 3]), 'rows')
unique(t_corners(:,[1 3]), 'rows')
unique(t_corners(:,[1 2]), 'rows')


corners_   = [corners; corners(1,:)];
t_corners_ = [t_corners; t_corners(1,:)];

figure(1)
plot3(  corners_(:,1),   corners_(:,2),   corners_(:,3)), grid on, hold on
plot3(t_corners_(:,1), t_corners_(:,2), t_corners_(:,3))
xlabel('x-Axis (m)'), ylabel('y-Axis (m)'), zlabel('z-Axis (m)')
arrow_length = 0.4;
quiver3(0, 0, 0, arrow_length, 0, 0, 'LineWidth', 2, 'AutoScale', 'off', 'color', [1 0 0])
quiver3(0, 0, 0, 0, arrow_length, 0, 'LineWidth', 2, 'AutoScale', 'off', 'color', [0 0.5 0])
quiver3(0, 0, 0, 0, 0, arrow_length, 'LineWidth', 2, 'AutoScale', 'off', 'color', [0 0 1]), hold off
axis equal