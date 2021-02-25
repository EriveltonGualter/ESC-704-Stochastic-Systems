close all
clear all;
clc;
       
addpath('acrobot'); addpath('Data'); addpath('helpers');

load('Data/results4.mat'); close all; p = param; p.iteration = 500;

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

disp('Press Enter to Continue ...');
pause();
clc

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

disp('Press Enter to Continue ...');
pause();
clc

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

disp('Press Enter to Continue ...');
pause();
clc

%% Simulate Acrobot
simAcrobot(X_sys,p)

disp('Press Enter to Continue ...');
pause();
clc

%% Simulate Acrobot (Show Opt Traject)
plotAcrobotOptTraj

disp('Press Enter to Continue ...');
pause();
clc

%% Load Pendulum Problem
clear all

load('Data/results_pendulum1.mat'); close all;

disp('Press Enter to Continue ...');
pause();
clc

cplot(x2,y2,'-','linewidth',3)


disp('Completed!');



