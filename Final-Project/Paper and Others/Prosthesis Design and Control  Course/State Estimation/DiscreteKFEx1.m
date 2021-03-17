function DiscreteKFEx1(N)

% Optimal State Estimation, by Dan Simon
%
% Discrete time Kalman filter for position estimation of a Newtonian system.
% This example illustrates the effectiveness of the Kalman filter for state
% estimation. It also shows how the variance of the estimation error 
% propagates between time steps and decreases as each measurement is processed.
% INPUTS: N = number of time steps.

if ~exist('N', 'var')
    N = 6;
end

T = 5; % time between measurements
sigma = 30; % position measurement standard deviation
R = sigma^2;
P0plus = diag([100, 10, 1]); % initial state estimate uncertainty

Q = diag([0, 0, 0]);
Q = diag([0, 0, 1]);

H = [1 0 0];
F = [1 T T*T/2; 0 1 T; 0 0 1]; % state transition matrix
x = [1; 1; 1]; % initial state
xhatplus = x; % initial state estimate

posArray = zeros(N+1, 1);
posArray(1) = x(1);
posEstArray = zeros(N+1, 1);
posEstArray(1) = xhatplus(1);
zArray = zeros(N+1, 1);
Pplus = P0plus;
Varminus = zeros(N+1, 1);
Varplus = zeros(N+1, 1);
Varplus(1) = P0plus(1,1);

for k = 1 : N
    % Simulate the system and measurement
    x = F * x;
    z = H * x + sigma * randn;
    % Estimate the state
    Pminus = F * Pplus * F' + Q;
    K = Pminus * H' / (H * Pminus * H' + R);
    xhatminus = F * xhatplus;
    xhatplus = xhatminus + K * (z - H * xhatminus);
    Pplus = (eye(3) - K * H) * Pminus * (eye(3) - K * H)' + K * R * K';
    % Save data for plotting
    posArray(k+1) = x(1);
    posEstArray(k+1) = xhatplus(1);
    zArray(k+1) = z;
    Varminus(k+1) = Pminus(1,1);
    Varplus(k+1) = Pplus(1,1);
end

% Plot the results
close all;
k = 0 : N;
figure; hold on; grid
plot(k, zArray-posArray, 'r:', 'LineWidth', 2);
plot(k, posEstArray-posArray, 'b-', 'LineWidth', 2);
set(gca,'FontSize',12); set(gcf,'Color','White'); set(gca,'Box','on');
xlabel('time step'); ylabel('position');
legend('measurement error', 'estimation error');

figure; hold on; grid
for k = 1 : N
    plot([k-1 k], [Varplus(k) Varminus(k+1)], 'LineWidth', 2);
    plot([k k], [Varminus(k+1) Varplus(k+1)], 'LineWidth', 2);
end
set(gca,'FontSize',12); set(gcf,'Color','White'); set(gca,'Box','on');
xlabel('time step');
ylabel('position estimation error variance');