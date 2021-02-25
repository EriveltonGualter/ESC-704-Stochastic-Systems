function dz = acrobotDynamics(z,u)

% Parameters
m1 = 1;             % Mass, kg
m2 = 1;             % Mass, kg
l1 = 1;             % Lenght, m
l2 = 1;             % Lenght, m
r1 = 0.5;           % Distance CG1, m
r2 = 0.5;           % Distance CG2, m
I1 = m1*l1^2/12;    % Moment of Inertia 1, m^3
I2 = m2*l2^2/12;    % Moment of Inertia 1, m^3
g  = 9.81;

dz = compute_accel(I1,I2,z(1),z(3),z(2),z(4),u(1),u(2),g,l1,m1,m2,r1,r2);