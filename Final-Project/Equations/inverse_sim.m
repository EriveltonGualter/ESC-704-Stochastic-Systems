clear all
close all
clc


dt = 1e-3;
t = 0:dt:4;

amp_ap = -40;
amp_bt =  40;
vx = 1;
w = 2*pi/t(end);

apb = amp_ap*cos(w*t)-20;
btb = amp_bt*cos(w*t)+60;
xb  = vx*t;

apbdot = -amp_ap*w*sin(w*t);
btbdot = -amp_bt*w*sin(w*t);
xbdot  = vx*ones(size(t));

apbddot = -amp_ap*w^2*cos(w*t);
btbddot = -amp_bt*w^2*cos(w*t);
xbddot  = zeros(size(t));

% figure;
% subplot(311); plot(t,xb);
% subplot(312); plot(t,apb);
% subplot(313); plot(t,btb);


qd = [xb; btb; apb; xbdot; btbdot; apbdot; xbddot; btbddot; apbddot];

q0 = qd(1:6,1);

% Gains
w = 2;
kp = diag([0 w^2 w^2]);
kd = diag([0 2*w 2*w]);

% 
% fun = @(t,q) sim_propulsion(t,q,qb,K);
% 
% sol = ode45(fun, t, q0);

q = q0;
for k=1:length(t)
    u(:,k) = inv_dyn_control([],q(:,k),qd(:,k),kp,kd);
    q(:,k+1) = q(:,k) + dynamicRK4(dt,q(:,k),u(:,k))*dt;
end

figure;
for i=1:6
    if i<4
        subplot(6,1,i); hold on; plot(q(i,:)); plot(qd(i,:)); ylim([min(qd(i,:)) max(qd(i,:))])
    else
        subplot(6,1,i); plot(u(i-3,:))        
    end
%     if i<7
%         subplot(9,1,i); plot(q(i,:))
%     else
%         subplot(9,1,i); plot(u(i-6,:))
%     end
end

function x_new=dynamicRK4(delta,x,u)
    %use Ruku4 for discretization
    k1=ode_dyn(x,u);
    k2=ode_dyn(x+delta/2*k1,u);
    k3=ode_dyn(x+delta/2*k2,u);
    k4=ode_dyn(x+delta*k3,u);
    x_new=x+delta/6*(k1+2*k2+2*k3+k4);
end


function qdot = sim_propulsion(t,q,qd,K)
    u = inv_dyn_control(t,q,qd,K);
    
    qdot = ode_dyn(q,u);
end