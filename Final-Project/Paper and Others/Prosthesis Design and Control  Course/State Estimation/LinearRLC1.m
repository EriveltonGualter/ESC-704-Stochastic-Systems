function LinearRLC1

R = 1;
L = 1;
C = 1;
A = [-2/R/C, 1/C ; -1/L, 0];
B = [1/R/C ; 1/L];
C = [-1, 0];
D = 1;
dt = 0.01;

t = 0 : dt : 5;
N = length(t);
yArray = zeros(1, N);
x = [0 ; 0];
for i = 1 : N-1
    u = 1;
    xdot = A * x + B * u;
    x = x + xdot * dt;
    y = C * x + D * u;
    yArray(i+1) = y;
end

close all
plot(t, yArray)
xlabel('time (s)')
ylabel('output (V)')