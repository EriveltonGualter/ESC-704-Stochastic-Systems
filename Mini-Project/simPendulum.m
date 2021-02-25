clearvars -global
close all

% Parameters
l1 = 1;     % Lenght, m
l2 = 1;     % Lenght, m

% Simulation parameters
dt          = 0.01; %p.dt;              % Integration step size [s]
iteration   = 500;  %p.iteration;
t           = 0:dt:dt*iteration; % array time

p.m1 = 1;   % (kg) Cart mass
p.m2 = 1;   % (kg) pole mass
p.l  = 1;   % (m) pendulum (pole) length

% Convert states to cartesian positions:
pos = cartPolePosition(X_sys,p);
    
x1 = pos(1,:);
y1 = pos(2,:);
x2 = pos(3,:);
y2 = pos(4,:);

% Plotting parameters:
p.w = 0.6*p.l;  % Width of the cart
p.h = 0.4*p.l;  % Height of the cart
p.r = 0.1*p.l;  % Radius of the pendulum bob
p.rw = 0.1*p.l; % Radius of the wheel

% Compute the extents of the drawing, keeping everything in view
padding = 0.2*p.l;  %Free space around edges
xLow = min(min(x1 - 0.5*p.w,  x2 - p.r)) - padding;
xUpp = max(max(x1 + 0.5*p.w,  x2 + p.r)) + padding;
yLow = min(min(y1 - 0.5*p.h,  y2 - p.r)) - padding;
yUpp = max(max(y1 + 0.5*p.w,  y2 + p.r)) + padding;
%     extents = [xLow,xUpp,yLow,yUpp];
extents = [-3 3 -2 2];

% Create and clear a figure:
fig = figure; fig.Renderer = 'painters'; cla; hold on;  box on;
set(gcf,'DoubleBuffer','on');   

% Compute the verticies of a star, just for fun;
star = getStarVerticies(7,0.5);     % 7 verticies, spoke ratio of 0.5
p.star = 0.6*p.r*star;              % Rescale the star;

time = 0;
tic;
i = 1;
while time < t(end)
    % Compute the position of the system at the current real world time
    posDraw = interp1(t',pos',time')';

    % Redraw the image
    drawCartPole(time,posDraw,extents,p);

    % Update current time
    time = toc;
end 

%% Plots
% figure;
% plot(t,z,'LineWidth',2);
% legend('displacement','velocity','angle','angular velocity');
% xlabel('time (s)','Fontsize',14,'interpreter','latex');
% ylabel('States','Fontsize',14,'interpreter','latex');
% xlim([min(t) max(t)]);

