function RLS
% Recursive least squares example
N = 100; % number of time steps
H = [1 1 ; 1 -1]; % measurement matrix
R = eye(2); % measurement noise covariance
x = [1 ; 1]; % constant vector that we want to estimate
xhat = [1/2 ; 3/4]; % initial estimate
P = eye(2); % initial estimation error covariance
xhatArray = zeros(2, N);
xhatArray(:, 1) = xhat;
PArray = zeros(2, N);
PArray(:, 1) = diag(P);
for k = 1 : N-1
    y = H * x + sqrt(R) * randn(2, 1);
    K = P * H' / (H * P * H' + R); % estimator gain
    xhat = xhat + K * (y - H * xhat); % estimator update
    P = (eye(2) - K * H) * P * (eye(2) - K * H)' + K * R * K';
    xhatArray(:, k+1) = xhat;
    PArray(:, k+1) = diag(P);
end
close all

figure, hold on
plot(xhatArray(1,:), 'r-')
plot(xhatArray(2,:), 'b--')
xlabel('time step')
ylabel('estimate')
legend('x(1)', 'x(2)')

figure, hold on
plot(PArray(1,:), 'r-')
plot(PArray(2,:), 'b--')
xlabel('time step')
ylabel('estimation error variance')
legend('P(1,1)', 'P(2,2)')