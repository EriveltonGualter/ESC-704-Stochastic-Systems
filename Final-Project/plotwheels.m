function hw = plotwheels(RW,LW,t,time,i,hw)
hold on;

% subplot(211); hold on; plot(MK.RW.y); plot(MK.LW.y); 
% subplot(212); hold on; plot(MK.RW.z); plot(MK.LW.z);

if i == 0
    [px,py,pz] = getSegment(RW, LW,t,time,i);
    RW.xi = px(1);
    RW.yi = py(1);
    RW.zi = pz(1);
    LW.xi = px(2);
    LW.yi = py(2);
    LW.zi = pz(2);
    N1 = [RW.xi-LW.xi, RW.yi-LW.yi, RW.zi-LW.zi];
    n1 = N1.'; n1(3)=0;
else
    N1 = [RW.x-LW.x, RW.y-LW.y, RW.z-LW.z];
    n1 = N1(i,:).'; n1(3)=0;
    
    RW.xi = RW.x(i);
    RW.yi = RW.y(i);
    LW.xi = LW.x(i);
    LW.yi = LW.y(i);
end

r = 299;

% Original points, original plane
t = linspace(0,2*pi);
x = r*cos(t);
y = r*sin(t);
z = 0*t;
pnts = [x;y;z];

% unit normal for original plane
n0 = [0;0;1]; 
n0 = n0/norm(n0);

% unit normal for plane to rotate into 
% plane is orthogonal to n1... given by equation
% n1(1)*x + n1(2)*y + n1(3)*z = 0
% n1 = [1;1;1]; 
n1 = n1/norm(n1); 

% theta is the angle between normals
c = dot(n0,n1) / ( norm(n0)*norm(n1) ); % cos(theta)
s = sqrt(1-c*c);                        % sin(theta)
u = cross(n0,n1) / ( norm(n0)*norm(n1) ); % rotation axis...
u = u/norm(u); % ... as unit vector

C = 1-c;

% the rotation matrix
R = [u(1)^2*C+c, u(1)*u(2)*C-u(3)*s, u(1)*u(3)*C+u(2)*s
    u(2)*u(1)*C+u(3)*s, u(2)^2*C+c, u(2)*u(3)*C-u(1)*s
    u(3)*u(1)*C-u(2)*s, u(3)*u(2)*C+u(1)*s, u(3)^2*C+c];

% Rotated points
newPnts = R*pnts;

% newPnts = newPnts - N1(1,:).';

if isempty(hw)
    hw(1) = plot3(newPnts(1,:)+RW.xi,newPnts(2,:)+RW.yi,newPnts(3,:)+r,'linewidth',2);
    hw(2) = plot3(newPnts(1,:)+LW.xi,newPnts(2,:)+LW.yi,newPnts(3,:)+r,'linewidth',2);
else
    set(hw(1),'xData',newPnts(1,:)+RW.xi,'yData',newPnts(2,:)+RW.yi,'zData',newPnts(3,:)+r);
    set(hw(2),'xData',newPnts(1,:)+LW.xi,'yData',newPnts(2,:)+LW.yi,'zData',newPnts(3,:)+r);
end

    
    
    
    
    