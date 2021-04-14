% clear all
% close all

% REJ2(:,1)); ylabel('Angle'); title('Alpha Left'); 
%     subplot(524); hold on; box on; plot(tswl, -LEJ2(:,1)); ylabel('Angle'); title('Alpha Right'); 
%     subplot(525); hold on; box on; plot(tswr, 1.8*RSJ); ylabel('Angle'); title('Beta Left')
%     subplot(526); hold on; box on; plot(tswl, 1.8*LSJ)

x = x_wr;
dx = dx_wr;
alph = 1.8*RSJ.'; 
beta = -LEJ2(:,1); 
dalph = diff([alph; alph(end)])/dt;
dbeta = diff([beta; beta(end)])/dt;

qref(:,1) = x;
qref(:,2) = alph;
qref(:,3) = beta;
qref(:,4) = dx;
qref(:,5) = dalph;
qref(:,6) = dbeta;

uref(:,1) = swl2.Fx'; 
uref(:,2) = swl2.Fy';
uref(:,3) = swl2.Mx';
uref(:,4) = swl2.My';

delta = 1/f1;
%%
clc
% CasADi v3.5.5
import casadi.*

qmax = max(qref);
qmin = min(qref);
umax = max(uref);
umin = min(uref);

t = SX.sym('t');
q = SX.sym('q',1,6);
u = SX.sym('u',1,4);

nx = length(q);
nu = length(u);
np = 13; 

N = length(qref)-1;

% Decision variables (controls)
U = SX.sym('U',nu,N);       
% Size nx+N*(nx+nu) correspond to trajectory and control, remaining correspong to the weights
P = SX.sym('P',np); 
% States over the optimization problem.
X = SX.sym('X',nx,(N+1));   

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
    % Control
    con = U(:,k);

    st_next_RK4 = dynamicRK4(delta,st,con);
    g = [g;st_next-st_next_RK4];
    
    % Desired Control
%     con_ref = P(nx*(k+1)+nu*(k-1)+1:nx*(k+1)+nu*(k-1)+nu);
    % Desired Trajectory
    st_ref = P(nx*k+nu*(k-1)+1:nx*(k+1)+nu*(k-1));
    
    % calculate obj
    obj = obj + runningcosts(st, con, qref(k,:).', uref(k,:).');
end

% make the decision variables one column vector
% OPT_variables = [reshape(X,nx*(N+1),1);reshape(U,nu*N,1)];
OPT_variables = [reshape(X,nx*(N+1),1);reshape(U,nu*N,1)];

% nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

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

opts.ipopt.linear_solver = 'MA86';

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);
% https://groups.google.com/g/casadi-users/c/28p-BHjgXTc/m/cXyqi48NBQAJ

args = struct;

args.lbg(1:nx*(N+1)) = 0;  % -1e-20  % Equality constraints
args.ubg(1:nx*(N+1)) = 0;  % 1e-20   % Equality constraints

% inequality constraints (state constraints)
args.lbx(1:nx:nx*(N+1),1) = qmin(1); 
args.ubx(1:nx:nx*(N+1),1) = qmax(1); 
args.lbx(2:nx:nx*(N+1),1) = qmin(2);
args.ubx(2:nx:nx*(N+1),1) = qmax(2);
args.lbx(3:nx:nx*(N+1),1) = qmin(3);
args.ubx(3:nx:nx*(N+1),1) = qmax(3);
args.lbx(4:nx:nx*(N+1),1) = qmin(4);
args.ubx(4:nx:nx*(N+1),1) = qmax(4);
args.lbx(5:nx:nx*(N+1),1) = qmin(5);
args.ubx(5:nx:nx*(N+1),1) = qmax(5);
args.lbx(6:nx:nx*(N+1),1) = qmin(6);
args.ubx(6:nx:nx*(N+1),1) = qmax(6);

% control constraints
args.lbx(nx*(N+1)+1:nu:nx*(N+1)+nu*N,1) = umin(1);
args.ubx(nx*(N+1)+1:nu:nx*(N+1)+nu*N,1) = umax(1);
args.lbx(nx*(N+1)+2:nu:nx*(N+1)+nu*N,1) = umin(2);
args.ubx(nx*(N+1)+2:nu:nx*(N+1)+nu*N,1) = umax(2);
args.lbx(nx*(N+1)+3:nu:nx*(N+1)+nu*N,1) = umin(3);
args.ubx(nx*(N+1)+3:nu:nx*(N+1)+nu*N,1) = umax(3);
args.lbx(nx*(N+1)+4:nu:nx*(N+1)+nu*N,1) = umin(4);
args.ubx(nx*(N+1)+4:nu:nx*(N+1)+nu*N,1) = umax(4);

%% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
close all

fig1 = figure; fig1.Renderer = 'painters'; 
fig2 = figure; fig2.Renderer = 'painters'; 

% for is=1:1%length(solution)
import casadi.*
% close all
% Weights
Q = [1 0 .1 0];
R = [1 1];
kp = 50;
kd = 5;
lambds = [1 1 0 0 0];

% initial value of the optimization variables
x0 = [0; 0.4; 0; 0];    % initial condition.
u0 = zeros(N,nu);       % two control inputs 
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables
args.x0  = [reshape(X0',nx*(N+1),1); reshape(u0',nu*N,1)];

%%% Trajectory 
Amp = deg2rad(30); 
Bias = Amp; 
f = 1/3;
tol = 0:dt:(N-1)*dt;
zref = zeros(4,N);
zref(1,:) = Amp*sin(2*pi*f*tol) + Bias; 
zref(3,:) = Amp*cos(2*pi*f*tol)*2*pi/4; 

% Parsing Trajectory to Args
args.p(1:nx) = x0; 
for k = 1:N 
    % Control ref (TODO: remove)
    u1_ref = 0; 
    u2_ref = 0;
    u3_ref = 0;
    
    args.p(nx*k+nu*(k-1)+1:nx*(k+1)+nu*(k-1))  = zref(:,k).';
    args.p(nx*(k+1)+nu*(k-1)+1: nx*(k+1)+nu*k) = [u1_ref, u2_ref];
end

% Parse Q and R weights
args.p(nx+N*(nx+nu)+1: nx+N*(nx+nu)+nq)             = Q;
args.p(nx+N*(nx+nu)+nq+1: nx+N*(nx+nu)+nq+nr)       = R;

% Parse Lambdas
args.p(nx+N*(nx+nu)+nq+nr+1: nx+N*(nx+nu)+nq+nr+nl) = lambds;

% PD gains
args.p(nx+N*(nx+nu)+nq+nr+nl+1) = kp;
args.p(nx+N*(nx+nu)+nq+nr+nl+2) = kd;
    
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
    
% Get External Moment from PD
Mext = kp*(xol(1:end-1,1)-zref(1,:)') + kd*(xol(1:end-1,3)-zref(3,:)');
uplot = [uol Mext];

disp(['External Power: ', num2str(trapz(tol, Mext))])
disp(['rms: ',num2str(rms(xol(1:end-1,1)-zref(1,:)'))])

ystr = {'q','x','qd','xd'};
% if exist('fig1','var')
%     if ishandle(fig1)
%         set(0, 'CurrentFigure', fig1); 
%         clf;
%     else
%         fig1 = figure; fig1.Renderer = 'painters'; set(gcf,'color','w');
%     end
% else
%     fig1 = figure; fig1.Renderer = 'painters'; set(gcf,'color','w');
% end
set(0, 'CurrentFigure', fig1); clf;
for i=1:size(xol,2)
    subplot(size(xol,2),1,i); hold on; box on; plot(tol,xol(1:end-1,i));
    if i==1 || i==3
        subplot(4,1,i); hold on; box on; plot(tol,zref(i,:))
    end
    plot(tol,Z_max(i)*ones(size(tol)),'r')
    plot(tol,Z_min(i)*ones(size(tol)),'r')
    
    ylim([Z_min(i) Z_max(i)])
    xlim([min(tol) max(tol)])
    ylabel(ystr{i})
end

ystr = {'Taua','Taub','Mext'};
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
for i=1:size(uplot,2)
    subplot(size(uplot,2),1,i); hold on; box on; 
    plot(tol,uplot(:,i));
    
    if i<3
        plot(tol,U_max(i)*ones(size(tol)),'r')
        plot(tol,U_min(i)*ones(size(tol)),'r')
    end
    
    xlim([min(tol) max(tol)]);
    ylabel(ystr{i});
end

% [dVs, Pe, e] = energyBalance(tol,xol,uol,zref,param);
% disp(dVs)
disp(['Delta Voltage: ', num2str(energyBalance(tol,xol(1:end-1,:),uol,zref,param,0))])
disp('----------------------------')
% end

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

function x_new=dynamicRK4(delta,x,u)
    %use Ruku4 for discretization   
    k1=ode_dyn(0,x,u);
    k2=ode_dyn(0,x+delta/2*k1,u);
    k3=ode_dyn(0,x+delta/2*k2,u);
    k4=ode_dyn(0,x+delta*k3,u);
    x_new=x+delta/6*(k1+2*k2+2*k3+k4);
end

function cost = runningcosts(x, u, xref, uref)
    Q = 1;
    R = 1;
    cost = (x-xref)'*Q*(x-xref) + (u-uref)'*R*(u-uref); 
end
