clear all
close all
clc


load('DT1.mat');

for i=1:length(uol)
   xoldot(:,i) = ode_dyn(xol(i,:),uol(i,:));
end

stry = {'x','bt','ap','xdot','btdot','apdot','xddot','btddot','apddot'};
for i=1:9
    if i <= 6
        subplot(9,1,i); hold on; plot(tol, xol(1:end-1,i));
    end
    if i>3
        subplot(9,1,i); hold on; plot(tol, xoldot(i-3,:));
    end
    ylabel(stry{i})
end

