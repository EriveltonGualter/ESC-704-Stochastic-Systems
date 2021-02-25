function drawCartPole(time,pos,extents,p)
% drawCartPole(time,pos,extents,p)
%
% INPUTS:
%   time = [scalar] = current time in the simulation
%   pos = [4, 1] = [x1;y1;x2;y2]; = position of [cart; pole]
%   extents = [xLow, xUpp, yLow, yUpp] = boundary of the draw window
%
% OUTPUTS:
%   --> A pretty drawing of a cart-pole
%

% Make plot objects global so that they don't get erased
global florHandle cartHandle poleHandle bobHandle starHandle ... 
    r_wheelHandle l_wheelHandle;

cartColor = [0.2, 0.7, 0.2];   % [R, G, B]
poleColor = [0.3, 0.2, 0.7];

%%%% Unpack the positions:
x1 = pos(1);
y1 = pos(2);
x2 = pos(3);
y2 = pos(4);

% Title and simulation time:
title(sprintf('t = %2.2f% s',time));


% Draw the rail that the cart-pole travels on
if isempty(florHandle)
    florHandle = plot(extents(1:2),[-p.h/2-2*p.rw,-p.h/2-2*p.rw],'k-','LineWidth',2);
end

% p.w = 0.6*p.l;  %Width of the cart
% p.h = 0.4*p.l;  %Height of the cart
% p.r = 0.1*p.l;  % Radius of the pendulum bob
% p.rw = 0.1*p.l; % Radius of the wheel



% Draw the cart:
if isempty(cartHandle)
    cartHandle =  rectangle(...
        'Position',[x1-0.5*p.rw, y1-0.5*p.rw, p.rw, p.rw],...
        'Curvature',0.5*[p.rw,p.rw]/p.l,...
        'LineWidth',2,...
        'FaceColor',cartColor,...
        'EdgeColor',0.5*cartColor);  %Darker version of color
else
    set(cartHandle,...
        'Position',[x1-0.5*p.w, y1-0.5*p.h, p.w, p.h]);
end

% Draw the rigth wheel:
if isempty(r_wheelHandle)
    r_wheelHandle = rectangle(...
        'Position',[x1+p.w/2-p.rw*2, y1-p.h/2-p.rw*2, 2*p.rw, 2*p.rw],...
        'Curvature',[1,1],...   % <-- Draws a circle...
        'LineWidth',2,...
        'FaceColor',poleColor,...
        'EdgeColor',0.5*poleColor);  %Darker version of color
else
    set(r_wheelHandle,...
        'Position',[x1+p.w/2-p.rw*2, y1-p.h/2-p.rw*2, 2*p.rw, 2*p.rw]);
end

% Draw the left wheel:
if isempty(l_wheelHandle)
    l_wheelHandle = rectangle(...
        'Position',[x1-p.w/2, y1-p.h/2-p.rw*2, 2*p.rw, 2*p.rw],...
        'Curvature',[1,1],...   % <-- Draws a circle...
        'LineWidth',2,...
        'FaceColor',poleColor,...
        'EdgeColor',0.5*poleColor);  %Darker version of color
else
    set(l_wheelHandle,...
        'Position',[x1-p.w/2, y1-p.h/2-p.rw*2, 2*p.rw, 2*p.rw]);
end

% Draw the pole:
if isempty(poleHandle)
    poleHandle = plot([x1,x2], [y1, y2],...
        'LineWidth',4,...
        'Color',0.8*poleColor);
else
    set(poleHandle,...
        'xData',[x1,x2],...
        'yData',[y1,y2]);
end


% Draw the bob of the pendulum:
if isempty(bobHandle)
    bobHandle = rectangle(...
        'Position',[x2-p.r, y2-p.r, 2*p.r, 2*p.r],...
        'Curvature',[1,1],...   % <-- Draws a circle...
        'LineWidth',2,...
        'FaceColor',poleColor,...
        'EdgeColor',0.5*poleColor);  %Darker version of color
else
    set(bobHandle,...
        'Position',[x2-p.r, y2-p.r, 2*p.r, 2*p.r]);
end


% Draw a star for fun:
if isempty(starHandle)
    starHandle = patch(...
        x1 + p.star(1,:),...
        y1 + p.star(2,:),...
        'r');
else
    set(starHandle,...
        'XData',x1 + p.star(1,:),...
        'YData',y1 + p.star(2,:));
end


% Format the axis so things look right:
axis equal; axis(extents); axis on;      %  <-- Order is important here

% Push the draw commands through the plot buffer
drawnow;

end