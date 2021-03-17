function RegenerationKalman

k = 20;
JG = 0.088e-4;
n = 70;
JM = 1.29e-4;
JT = JG + n^2 * JM;
a = 0.054;
R = 0.0821;
L = 0.001;
u = 68;
C = 200;

A = [0, 1, 0, 0 ; -k/JT, 0, -n*a/JT, 0 ; 0, a*n/L, -R/L, -u/L ; 0, 0, u/C, 0];
B = [0 ; 1/JT ; 0 ; 0];
Q = diag([0, 0.01, 0.01, 0.01]);
C = [1, 0, 0, 0];
R = 0.01; % Try different values of R and see what happens!

load 'WinterData.mat'; % Load t (time) and M (input torque)
dt = 0.0001; % Try a larger step size and see what happens!
tfine = t(1) : dt : t(end); % Get better resolution for the input data
Mfine = spline(t, M, tfine);

N = length(tfine);
xArray = zeros(4, N);
x = [0 ; 0 ; 0 ; 0];
xArray(:, 1) = x;
P = 0.01 * eye(4);
xhat = x + sqrt(P) * randn(4, 1);
xhatArray(:, 1) = xhat;
KArray = zeros(4, N-1);
for i = 1 : N-1
    xdot = A * x + B * Mfine(i) + sqrt(Q) * randn(4, 1);
    x = x + xdot * dt;
    y = C * x + sqrt(R) * randn;
    xArray(:, i+1) = x;
    % Kalman filter
    K = P * C' / R;
    xhatdot = A * xhat + B * Mfine(i) + K * (y - C * xhat);
    xhat = xhat + xhatdot * dt;
    Pdot = -P * C' / R * C * P + A * P + P * A' + Q;
    P = P + Pdot * dt;
    xhatArray(:, i+1) = xhat;
    KArray(:, i) = K;
end

close all

figure
plot(tfine, xArray(1,:), 'r-', tfine, xhatArray(1,:), 'b--')
xlabel('time (s)')
ylabel('knee angle (rad)')
legend('true', 'estimate')
RMSError = sqrt(sum((xArray(1,:)-xhatArray(1,:)).^2)/N);
disp(['RMS knee angle estimation error = ', num2str(RMSError), ' rad'])

figure
plot(tfine, xArray(2,:), 'r-', tfine, xhatArray(2,:), 'b--')
xlabel('time (s)')
ylabel('knee velocity (rad/sec)')
legend('true', 'estimate')
RMSError = sqrt(sum((xArray(2,:)-xhatArray(2,:)).^2)/N);
disp(['RMS knee velocity estimation error = ', num2str(RMSError), ' rad/sec'])

figure
plot(tfine, xArray(3,:), 'r-', tfine, xhatArray(3,:), 'b--')
xlabel('time (s)')
ylabel('motor current (Amps)')
legend('true', 'estimate')
RMSError = sqrt(sum((xArray(3,:)-xhatArray(3,:)).^2)/N);
disp(['RMS motor current estimation error = ', num2str(RMSError), ' A'])

figure
plot(tfine, xArray(4,:), 'r-', tfine, xhatArray(4,:), 'b--')
xlabel('time (s)')
ylabel('capacitor voltage (V)')
legend('true', 'estimate')
RMSError = sqrt(sum((xArray(4,:)-xhatArray(4,:)).^2)/N);
disp(['RMS capacitor voltage estimation error = ', num2str(RMSError), ' V'])

figure
plot(tfine(2:end),KArray(1,:), tfine(2:end),KArray(2,:), tfine(2:end),KArray(3,:), tfine(2:end),KArray(4,:))
xlabel('time (s)')
ylabel('Kalman gain')
legend('K(1)', 'K(2)', 'K(3)', 'K(4)')
