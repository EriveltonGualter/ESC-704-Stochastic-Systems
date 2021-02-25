% Adapted code from stackoverflow
% https://stackoverflow.com/questions/14399689/matplotlib-drawing-lines-between-points-ignoring-missing-data

% Parameters
l1 = 1;     % Lenght, m
l2 = 1;     % Lenght, m

% Simulation parameters
% dt          = p.dt;              % Integration step size [s]
% iteration   = p.iteration;
dt = 0.01;
iteration = 500;
t           = 0:dt:dt*iteration; % array time

% Unpack Angles
q1 = X_sys(1,:);
q2 = X_sys(3,:);

% Forward Kinematics
x1 = l1*cos(q1);
y1 = l1*sin(q1);
x2 = l1*cos(q1) + l2*cos(q1+q2);
y2 = l1*sin(q1) + l2*sin(q1+q2);

pos = [x1; y1; x2; y2];

skp = 10;
x = [zeros(size(x1)); x1; x2]; 
y = [zeros(size(y1)); y1; y2]; 

x = x(:,1:skp:end);
y = y(:,1:skp:end);

% "Buffer" size, number of historic lines to keep, and governs the 
% corresponding fade increments.
n = size(x,2);
nFade = 100;

% Set up some demo values for plotting around a circle
dt = 0.05; 
a = 0:dt:2*pi+(dt*nFade); 
c = colormap(parula(size(x,2)));

% Initialise the figure, set up axes etc
extents = [-2.2 2.2 -2.2 2.2];

f = figure(1); clf; axis equal; axis(extents); axis on; box on;
f.Renderer = 'painters';

% Draw all of the lines, initially not showing because NaN vs NaN
lines = arrayfun( @(x)line(NaN,NaN), 1:nFade, 'uni', 0 );
% Set up shorthand for recolouring all the lines
recolour = @(lines) arrayfun( @(x) set( lines{x},'Color',ones(1,3)*(1-x/nFade) ), 1:nFade );
% recolour = @(lines) arrayfun( @(x) set( lines{x},'Color',c(x,:)), 1:nFade );

j = 1;
for ii = 1:n
    % Shift the lines around so newest is at the start
    lines = [ lines(end), lines(1:end-1) ]; 
    % Overwrite x/y data for oldest line to be newest line
    set( lines{1}, 'XData', x(:,ii), 'YData', y(:,ii) );
    % Update all colours
    recolour( lines );           
    
    % Pause for animation           
    pause(0.01);   
end

hold on;

for ii = 1:n
    if any([1 5 10 15 22 30 40 51] == ii)
        hold on
        clearvars -global
        draw(x(2,ii), y(2,ii), x(3,ii), y(3,ii), c(ii,:));
    end
    
    % Pause for animation           
    pause(0.01);
end

hold on;

x = [zeros(size(x1)); x1; x2]; x = x(3,:).';
y = [zeros(size(y1)); y1; y2]; y = y(3,:).';

cplot(x,y,'-','linewidth',3)

