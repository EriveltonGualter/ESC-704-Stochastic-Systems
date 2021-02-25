close all
clear all;
clc;

% Implementation inspired from paper:
%   Entropy code from : Williams, G., Drews, P., Goldfain, B., Rehg, J.
%   M., & Theodorou, E. A. (2016, May). Aggressive driving with model
%   predictive path integral control. In 2016 IEEE International
%   Conference on Robotics and Automation (ICRA) (pp. 1433-1440). IEEE.
        
addpath('acrobot'); addpath('Data'); addpath('helpers');

% MPC Parameters
K           = 1000;
N           = 100;
p.iteration = 500;
p.dt        = 0.01;

% Variance and Lamda
p.lambda    = 250;
p.variance  = 10;
p.R         = 1;

% Initiazation of variables
X_sys       = zeros(4,p.iteration+1);
U_sys       = zeros(2,p.iteration);
cost        = zeros(1,p.iteration);
cost_avg    = zeros(1,p.iteration);
x           = zeros(4,N);
u           = zeros(2,N);
delta_u     = zeros(N,K);
u_init      = 1;

% Initial State
x0          = [0 0 0 0];
X_sys(:,1)  = x0;

% Final State
p.xT        = [pi/2 0 0 0];

% Positive-Definite Matrix R
p.R         = diag([1e3 1e-3 10 0.1]);

for j = 1: p.iteration
    % Initialization of Cost for K Samples
    Stk = zeros(1,K);
    
    % Calculating cost for K samples and N finite horizone
    for k = 1:K
        x(:,1) = x0;
        for i = 1:N-1
            delta_u(i,k)= p.variance*(randn(1));
            x(:,i+1)    = x(:,i) + sode_function(x(:,i), u(:,i)+delta_u(i,k), p)*p.dt;
            Stk(k)      = Stk(k) + cost_function(x(:,i+1), u(:,i)+ delta_u(i,k),p);
        end
        delta_u(N,k) = p.variance*(randn(1));
    end
    
    % Average cost over p.iteration
    cost_avg(j) = sum(Stk)/K;
    
    % Updating the control input according to the expectency over K sample trajectories
    for i = 1:N
        u(:,i) = u(:,i) + totalEntropy(Stk(:), delta_u(i,:),p);
    end
    
    % Input to the system
    U_sys(:,j) = u(:,1);
    
    % Discretized plant 
    X_sys(:,j+1) = X_sys(:,j) + sode_function(X_sys(:, j), u(:,1), p)*p.dt;
    
    % Find Cost function
    cost(j+1) = cost_function(X_sys(:,j+1),(u(:,i)+delta_u(i,k)),p);
    
    % Updating the input
    for i = 1:N-1
        u(:,i) = u(:,i+1);
    end
    u(N) = u_init;
    
    % Expectency over trajectory in next step
    x0 = X_sys(:,j+1);

    if sum(sum(isnan(X_sys)))
        break
    end
end

%% Create State Response Figure
fig = figure; hold on; box on;
fig.Renderer = 'painters';

% Setting subplot titles
str_titles = ...
    {'Angle Position - Link 1', ...
    'Angular Velocity - Link 1', ...
    'Angle Position - Link 2', ...
    'Angular Velocity - Link 2'};

% Create subplot
for i=1:4
    subplot(4,1,i); plot(X_sys(i,:),'LineWidth',2);  xlim([0 500])
    title(str_titles(i),    'Interpreter','Latex', 'FontSize',12);
    xlabel('Interaction',   'Interpreter','Latex', 'FontSize',12);
end

% Save figure
% saveas(gcf,'states','epsc')

%% Create Control Figuere
fig = figure; hold on; box on;
fig.Renderer = 'painters';

% Setting subplot titles
str_titles = ...
    {'Torque control applied to joint 1', ...
    'Torque control applied to joint 2'};

% Create subplot
for i=1:2
    subplot(2,1,i); plot(U_sys(i,:),'LineWidth',2);  xlim([0 500])
    title(str_titles(i),    'Interpreter','Latex', 'FontSize',12);
    xlabel('Interaction',   'Interpreter','Latex', 'FontSize',12);
end

% Save figure
% saveas(gcf,'control','epsc')

%% Create pdf figure
fig = figure; hold on; box on;
fig.Renderer = 'painters';

% Reshape random pertubation 
du      = reshape(delta_u,1,numel(delta_u));

% Find Kernel smoothing function estimate for univariate and bivariate data
[f,xi]  = ksdensity(du);

% Create normal distribution distribution object 
pd          = makedist('Normal','mu',0,'sigma',10);

% Solve Probability density function
x           = linspace(-50,50,1e3);
pdf_normal  = pdf(pd,x);

% Plots: Histogram, Kernel smoothing, and pdf
histogram(du,'Normalization','probability');
plot(xi,f,'LineWidth',3); 
plot(x,pdf_normal,'LineWidth',2); 

% Plot seetings
title('Probability density function',   'Interpreter','Latex', 'FontSize',12);
xlabel('Control pertubation',           'Interpreter','Latex', 'FontSize',12);
legend('Histogram from RV', 'Kernel smoothing', 'pdf', 'Interpreter','Latex', 'FontSize',12);

% Save figure
% saveas(gcf,'pdf','epsc')

%% Simulate Acrobot
simAcrobot(X_sys,p)

%% Simulate Acrobot (Show Opt Traject)
plotAcrobotOptTraj


