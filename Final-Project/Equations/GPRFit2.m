close all
clear
clc

Script4

load('DT1.mat');
tol = swl2.t';

fig = figure; fig.Renderer = 'painters';
fig.Position(3) = fig.Position(3)*2;
fig.Position(4) = fig.Position(4)*2;

tiledlayout(2,2,'TileSpacing','compact')

xol = [x_wl -LEJ2(:,1) 1.8*LSJ'];
xold = diff(xol,2);
xol = horzcat(xol, [xold; xold(end-1:end,:)]);
uol = [swl2.Fx' swl2.Fy' swl2.Mx' swl2.My'];

for k=1:4
x_observed = xol(1:end,:);
% x_observed = xol(1:end-1,:);
y_observed1 = uol(:,k);
y_observed2 = y_observed1 + 1*randn(size(y_observed1));

gprMdl = fitrgp(x_observed,y_observed2);

[ypred,~,yint] = predict(gprMdl,x_observed);

nexttile
scatter(tol,y_observed1,'.')                 % Observed data points
hold on; box on;scatter(tol,y_observed2,'xr') % Observed data points
plot(tol,ypred,'g', 'linewidth',1)           % GPR predictions
patch([tol;flipud(tol)],[yint(:,1);flipud(yint(:,2))],'k','FaceAlpha',0.1); % Prediction intervals
hold off
title('GPR Fit of Noisy Observations')
legend({'Noise-free observations','Noisy observations','GPR predictions','95% prediction intervals'},'Location','best')
ylim([min(yint(:,1)) max(yint(:,2))]);
end

