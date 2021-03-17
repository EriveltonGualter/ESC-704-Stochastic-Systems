function LinearRLC1Discrete

R = 1;
L = 1;
C = 1;
A = [-2/R/C, 1/C ; -1/L, 0];
B = [1/R/C ; 1/L];
C = [0, L];
dt = 0.01;

F = expm(A*dt);
G = (F - eye(2)) / A * B;

t = 0 : dt : 5;
N = length(t);
yArray = zeros(1, N);
x = [0 ; 0];
for i = 1 : N-1
    u = 1;
    x = F * x + G * u;
    y = C * x;
    yArray(i+1) = y;
end

close all
plot(t, yArray)
xlabel('time (s)')
ylabel('output (V)')