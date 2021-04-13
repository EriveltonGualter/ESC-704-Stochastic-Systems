%NO REFERENCE CONTROL with PD
clear all
% close all
clc

% CasADi v3.5.5
import casadi.*

%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SETTING UP

% dt = 0.01; % sampling time [s]
dt = .01;
T = 4;
N = T/dt; % prediction horizon

z1_max = 10;    z1_min = 0;
z2_max = pi;    z2_min = -z2_max;
z3_max = pi;    z3_min = -z3_max;
z4_max = 10;    z4_min = -z4_max;
z5_max = 10;    z5_min = -z5_max;
z6_max = 10;    z6_min = -z6_max;
u1_max = 50;    u1_min = -u1_max;
u2_max = 50;	u2_min = -u2_max;
u3_max = 20;    u3_min = -u3_max;
u4_max = 20;    u4_min = -u4_max;

Z_min = [z1_min z2_min z3_min z4_min z5_min z6_min];
Z_max = [z1_max z2_max z3_max z4_max z5_max z6_max];
U_min = [u1_min u2_min u3_min u4_min];
U_max = [u1_max u2_max u3_max u4_max];

z1 = SX.sym('z1'); 
z2 = SX.sym('z2');
z3 = SX.sym('z3');
z4 = SX.sym('z4');
z5 = SX.sym('z5');
z6 = SX.sym('z6');
states = [z1;z2;z3;z4;z5;z6]; 
nx = length(states);

u1 = SX.sym('u1'); 
u2 = SX.sym('u2');
u3 = SX.sym('u3');
u4 = SX.sym('u4');
controls = [u1;u2;u3;u4];
nu = length(controls);

% length(Q) + length(R) = nx+nu 

% Number of Weights
nq = nx;    % Elements of Q
nr = nu;    % Elements of R
% Total input weights
nw = nq+nr;

% Decision variables (controls)
U = SX.sym('U',nu,N);       
% Size nx+N*(nx+nu) correspond to trajectory and control, remaining correspong to the weights
P = SX.sym('P',nx+N*(nx+nu)+nw); 
% States over the optimization problem.
X = SX.sym('X',nx,(N+1));   

dz = ode_dyn(states,controls);
FSX = Function('f',{states,controls},{dz});

%%% Initialization of variables
% Objective function
obj = 0;
% constraints vector
g = [];

%%% compute solution symbolically
% initial state
st  = X(:,1); 
% initial condition constraints
g = [g;st-P(1:nx)]; 
%%% Generate NLP
for k = 1:N
    % Current State
    st = X(:,k); 
    % Next State
    st_next = X(:,k+1); 
    % Desired Trajectory
    st_ref = P(nx*k+nu*(k-1)+1:nx*(k+1)+nu*(k-1));
    % Control
    con = U(:,k);
   
    % Weitghs
    Q = diag(P(nx+N*(nx+nu)+1: nx+N*(nx+nu)+nq));
    R = diag(P(nx+N*(nx+nu)+nq+1: nx+N*(nx+nu)+nq+nr));
               
    st_next_RK4 = dynamicRK4(FSX,dt,st,con);

    g = [g;st_next-st_next_RK4]; % compute constraints % new
    
    if k~=N
        con_next = U(:,k+1);
    else
        con_next = U(:,k);
    end
    % calculate obj
    obj = obj + runningcosts(st, con, con_next, st_ref, Q, R);
end

% make the decision variables one column vector
OPT_variables = [reshape(X,nx*(N+1),1);reshape(U,nu*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

%%
import casadi.*

opts = struct;
opts.ipopt.max_iter = 2000; %2000
opts.ipopt.print_level =3;%0,3
opts.print_time = 1;
opts.ipopt.acceptable_tol =1e-6;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

% opts.expand = true;
% opts.jit = true;
% opts.compiler = 'shell';
% opts.jit_options.flags = {'-O3'};
% opts.jit_options.verbose = true;

% opts.ipopt.linear_solver = 'MA86';

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);
% https://groups.google.com/g/casadi-users/c/28p-BHjgXTc/m/cXyqi48NBQAJ

args = struct;

args.lbg(1:nx*(N+1)) = 0;  % -1e-20  % Equality constraints
args.ubg(1:nx*(N+1)) = 0;  % 1e-20   % Equality constraints

% inequality constraints (state constraints)
args.lbx(1:nx:nx*(N+1),1) = z1_min; 
args.ubx(1:nx:nx*(N+1),1) = z1_max; 
args.lbx(2:nx:nx*(N+1),1) = z2_min; 
args.ubx(2:nx:nx*(N+1),1) = z2_max; 
args.lbx(3:nx:nx*(N+1),1) = z3_min;
args.ubx(3:nx:nx*(N+1),1) = z3_max;
args.lbx(4:nx:nx*(N+1),1) = z4_min; 
args.ubx(4:nx:nx*(N+1),1) = z4_max; 
args.lbx(5:nx:nx*(N+1),1) = z5_min; 
args.ubx(5:nx:nx*(N+1),1) = z5_max; 
args.lbx(6:nx:nx*(N+1),1) = z6_min; 
args.ubx(6:nx:nx*(N+1),1) = z6_max; 

% control constraints
args.lbx(nx*(N+1)+1:nu:nx*(N+1)+nu*N,1) = u1_min; 
args.ubx(nx*(N+1)+1:nu:nx*(N+1)+nu*N,1) = u1_max; 
args.lbx(nx*(N+1)+2:nu:nx*(N+1)+nu*N,1) = u2_min;
args.ubx(nx*(N+1)+2:nu:nx*(N+1)+nu*N,1) = u2_max;
args.lbx(nx*(N+1)+3:nu:nx*(N+1)+nu*N,1) = u3_min;
args.ubx(nx*(N+1)+3:nu:nx*(N+1)+nu*N,1) = u3_max;
args.lbx(nx*(N+1)+4:nu:nx*(N+1)+nu*N,1) = u4_min;
args.ubx(nx*(N+1)+4:nu:nx*(N+1)+nu*N,1) = u4_max;

% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------

fig1 = figure; fig1.Renderer = 'painters'; 
fig2 = figure; fig2.Renderer = 'painters'; 

% for is=1:1%length(solution)
import casadi.*
% close all
% Weights
Q = [1e2 1e4 1e4 0 0 0];
R = [1 1 1 1]*1e-3;

%%% Trajectory 
tol = 0:dt:(N-1)*dt;
amp_ap = -deg2rad(40);
amp_bt =  deg2rad(40);
vx = 1;
w = 2*pi;
apb = amp_ap*cos(w*tol)-deg2rad(20);
btb = amp_bt*cos(w*tol)+deg2rad(60);
xb  = vx*tol;
apbdot = -amp_ap*w*sin(w*tol);
btbdot = -amp_bt*w*sin(w*tol);
xbdot  = vx*ones(size(tol));
apbddot = -amp_ap*w^2*cos(w*tol);
btbddot = -amp_bt*w^2*cos(w*tol);
xbddot  = zeros(size(tol));
qd = [xb; btb; apb; xbdot; btbdot; apbdot; xbddot; btbddot; apbddot];
zref(1:6,:) = qd(1:6,:); 

% initial value of the optimization variables
x0 = zref(:,1);    % initial condition.
u0 = zeros(N,nu);       % two control inputs 
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables
args.x0  = [reshape(X0',nx*(N+1),1); reshape(u0',nu*N,1)];



% Parsing Trajectory to Args
args.p(1:nx) = x0; 
for k = 1:N 
    % Control ref (TODO: remove)
    u1_ref = 0; 
    u2_ref = 0;
    u3_ref = 0;
    u4_ref = 0;
    
    args.p(nx*k+nu*(k-1)+1:nx*(k+1)+nu*(k-1))  = zref(:,k).';
    args.p(nx*(k+1)+nu*(k-1)+1: nx*(k+1)+nu*k) = [u1_ref, u2_ref, u3_ref, u4_ref];
end

% Parse Q and R weights
args.p(nx+N*(nx+nu)+1: nx+N*(nx+nu)+nq)             = Q;
args.p(nx+N*(nx+nu)+nq+1: nx+N*(nx+nu)+nq+nr)       = R;
    
% Solver
tic
sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
                'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
toc

% Cost Function
% f = NLP_casadi_exo(lambds(1:2), zref.', solver, args, param,Q,R);

% Get controls only from the solution
uol = reshape(full(sol.x(nx*(N+1)+1:end))',nu,N)'; 

% Get OPTIMAL solution TRAJECTORY
xol= reshape(full(sol.x(1:nx*(N+1)))',nx,N+1)'; 
    

%%
set(0, 'CurrentFigure', fig1); clf;
for i=1:size(xol,2)
    subplot(size(xol,2),1,i); hold on; box on; plot(tol,xol(1:end-1,i));
    if i<4
        plot(tol,qd(i,:));
    end
%     if i==1 || i==3
%         subplot(4,1,i); hold on; box on; plot(tol,zref(i,:))
%     end
%     plot(tol,Z_max(i)*ones(size(tol)),'r')
%     plot(tol,Z_min(i)*ones(size(tol)),'r')
%     
%     ylim([Z_min(i) Z_max(i)])
%     xlim([min(tol) max(tol)])
%     ylabel(ystr{i})
end

% if exist('fig2','var')
%     if ishandle(fig2)
%         set(0, 'CurrentFigure', fig2); 
%         clf;
%     else
%         fig2 = figure; fig2.Renderer = 'painters'; set(gcf,'color','w');
%     end
% else
%     fig2 = figure; fig2.Renderer = 'painters'; set(gcf,'color','w');
% end
set(0, 'CurrentFigure', fig2); 
for i=1:size(uol,2)
    subplot(size(uol,2),1,i); hold on; box on; 
    plot(tol,uol(:,i));
    
%     if i<3
%         plot(tol,U_max(i)*ones(size(tol)),'r')
%         plot(tol,U_min(i)*ones(size(tol)),'r')
%     end
%     
%     xlim([min(tol) max(tol)]);
%     ylabel(ystr{i});
end

%% Assemply psolver
% psolver.args = args;
% psolver.param = param;
% psolver.nx = nx;
% psolver.nu = nu;
% psolver.nq = nq;
% psolver.nr = nr;
% psolver.nl = nl;
% psolver.N = N;
% psolver.Q = Q;
% psolver.R = R;
% psolver.kp = kp;
% psolver.kd = kd;
% psolver.lbds = [0 0];
% % 
% nvar = 3; lb = [0 0 0]; ub = [1 1 1];

% f = NLP_casadi_exo([1 1 1], zref, solver,psolver)

%% Functions

function x_new=dynamicRK4(f,delta,x,u)
    %use Ruku4 for discretization
    k1=f(x,u);
    k2=f(x+delta/2*k1,u);
    k3=f(x+delta/2*k2,u);
    k4=f(x+delta*k3,u);
    x_new=x+delta/6*(k1+2*k2+2*k3+k4);
end

function cost = runningcosts(x, u, u_nxt, x_eq, Q, R)
    % Provide the running cost   
    
    J_err = (x-x_eq)'*Q*(x-x_eq);
    J_udt = (u_nxt-u)'*R*(u_nxt-u);
    J_u   = u'*R*u;
    
    cost = J_err + J_udt + J_u;
end
