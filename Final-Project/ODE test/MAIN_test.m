clear all
close all
clc

dt = 1e-2;
t = 0:dt:10; 

A = 1;
w = 2*pi/4;
xref = A*sin(w*t);
dxref = A*w*cos(w*t);
ddxref = -A*w^2*sin(w*t);

kp = 50;
kd = 5;

XR = [xref' dxref' ddxref'];

x1(:,1) = [0 0];
for i=1:length(t)-1
    e = XR(i,1) - x1(1,i);
    ed = XR(i,2) - x1(2,i);
    u1(i) = kp*e + kd*ed;
    
    x1(:,i+1) = dynamicRK4(dt,x1(:,i),u1(i),1);
end

subplot(311); hold on; plot(t, x1(1,:)); plot(t, XR(:,1));
subplot(312); hold on; plot(t, x1(2,:)); plot(t, XR(:,2));
subplot(313); plot(t(1:end-1), u1)

%%
import casadi.*

x = SX.sym('x',2,1);
u = SX.sym('u');

nx = length(x);
nu = length(u);

N = length(u1);

Q = eye(2);
R = 1;

X = SX.sym('X',nx,(N+1));  
U = SX.sym('U',nu,N);     
P = SX.sym('P',1); 
obj = 0;
g = [];

st= X(:,1); 
con = U(:,1);
g = [g;st-x1(:,1); con-u1(1)]; 
for k=1:N
    st = X(:,k); 
    st_next = X(:,k+1); 
%     st_ref = P(nx*k+nu*(k-1)+1:nx*(k+1)+nu*(k-1));
    con = U(:,k);
    
    st_next_RK4 = dynamicRK4(dt,st,con,P);
    g = [g;st_next-st_next_RK4]; % compute constraints % new
    g = [g; x1(1,k)-st(1)];
    obj = obj + runningcosts(st, con, x1(:,k), u1(:,k), Q, R);
end

OPT_variables = [reshape(X,nx*(N+1),1);reshape(U,nu*N,1);P];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g);


opts = struct;
opts.ipopt.max_iter = 2000; %2000
opts.ipopt.print_level =3;%0,3
opts.print_time = 1;
opts.ipopt.acceptable_tol =1e-9;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);
%%

args.lbg = 0;  % -1e-20  % Equality constraints
args.ubg = 0;  % 1e-20   % Equality constraints

% inequality constraints (state constraints)
minx = min(x1');
maxx = max(x1');
maxu = max(u1');
minu = min(u1');
args.lbx(1:nx:nx*(N+1),1) = minx(1); 
args.ubx(1:nx:nx*(N+1),1) = maxx(1); 
args.lbx(2:nx:nx*(N+1),1) = minx(2); 
args.ubx(2:nx:nx*(N+1),1) = maxx(2); 
args.lbx(nx*(N+1)+1:nu:nx*(N+1)+nu*N,1) = minu(1); 
args.ubx(nx*(N+1)+1:nu:nx*(N+1)+nu*N,1) = maxu(1); 
args.lbx(nx*(N+1)+nu*N+1) = 0; 
args.ubx(nx*(N+1)+nu*N+1) = 1.1; 

% control constraints
% args.lbx(nx*(N+1)+1:nu:nx*(N+1)+nu*N,1) = u1_min; 
% args.ubx(nx*(N+1)+1:nu:nx*(N+1)+nu*N,1) = u1_max; 
% args.lbx(nx*(N+1)+2:nu:nx*(N+1)+nu*N,1) = u2_min;
% args.ubx(nx*(N+1)+2:nu:nx*(N+1)+nu*N,1) = u2_max;

% args.x0 = 1; 

import casadi.*

args.x0 = [reshape(x1,nx*(N+1),1);reshape(u1,nu*N,1);0];

sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx, 'lbg', args.lbg, 'ubg', args.lbg);

xol= reshape(full(sol.x(1:nx*(N+1)))',nx,N+1); 
uol = reshape(full(sol.x(nx*(N+1)+1:end-1))',nu,N);
pol = full(sol.x(end))


figure;
subplot(311); hold on; plot(t, xol(1,:)); plot(t, x1(1,:));
subplot(312); hold on; plot(t, xol(2,:)); plot(t, x1(2,:));
subplot(313); hold on; plot(t(1:end-1), uol); plot(t(1:end-1), u1); 

%% Functions

function x_new=dynamicRK4(delta,x,u,p)
    %use Ruku4 for discretization   
    k1=ode_test(0,x,u,p);
    k2=ode_test(0,x+delta/2*k1,u,p);
    k3=ode_test(0,x+delta/2*k2,u,p);
    k4=ode_test(0,x+delta*k3,u,p);
    x_new=x+delta/6*(k1+2*k2+2*k3+k4);
end

function cost = runningcosts(x, u, u_des, x_des, Q, R)
    cost = (x-x_des)'*Q*[1 0; 0 0]*(x-x_des);
    %+ (u-u_des)'*R*(u-u_des);
end
