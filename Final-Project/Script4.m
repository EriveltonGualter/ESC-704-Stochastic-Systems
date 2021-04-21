
addpath('helpers');

% trial = 2;
% if ~exist('swl2', 'var')
    run('Script2.m');
% end
dt = 1/f1;

% Assembly U
UL = [swl2.Fx' swl2.Fy' swl2.Fz' swl2.Mx' swl2.My' swl2.Mz' swl2.t'];
UR = [swr2.Fx' swr2.Fy' swr2.Fz' swr2.Mx' swr2.My' swr2.Mz' swr2.t'];

% [m,n] = size(UL);


% x = zeros(n-1,m+1);
% for k=1:m
%     x(:,k+1) = x(:,k) + ode_dyn(x(:,k),UL(k,1:end-1))*dt;
% end
% 
% for i=1:6
%     subplot(6,1,i); plot(0:dt:dt*m, x(i,:));
% end

x_wr = sqrt(sum([MK.RW.x MK.RW.y MK.RW.z].^2,2))*1e-3;
x_wl = sqrt(sum([MK.LW.x MK.LW.y MK.LW.z].^2,2))*1e-3;

dx_wr = diff([x_wr; x_wr(end)])/dt;
dx_wl = diff([x_wl; x_wl(end)])/dt;

% figure;
% subplot(221); plot(t, x_wr);
% subplot(222); plot(t, x_wl);
% subplot(223); plot(t, dx_wr);
% subplot(224); plot(t, dx_wl);

%%
figure;
clf
box on; hold on;
plot(t, -REJ2(:,1));
plot(t, 1.8*RSJ);
drawnow

%% ALL PLOTS
figure; 
plotDTforRegression(swl2, swr2, REJ2, LEJ2, RSJ, LSJ, x_wr, x_wl, dx_wr, dx_wl)

% saveas(gcf,['plotresults/data_trial',num2str(idfig),'.pdf'])

set(gcf,'renderer','Painters')
% saveas(gcf,['plotresults/gpr_trial',num2str(idfig),'.pdf'])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,['plotresults/data_trial',num2str(idfig)],'-dpdf','-r0')


clc
% add some noise
swlnoise    = addNoise(swl2, 5);
swrnoise    = addNoise(swr2, 5);
REJnoise    = addNoise(REJ2, 5);
LEJnoise    = addNoise(LEJ2, 5);
RSJnoise    = addNoise(RSJ, 5);
LSJnoise    = addNoise(LSJ , 5);
x_wrnoise   = addNoise(x_wr, 0.1);
x_wlnoise   = addNoise(x_wl, 0.1);
dx_wrnoise  = addNoise(dx_wr, 0.1);
dx_wlnoise  = addNoise(dx_wl, 0.1);

% figure; 
plotDTforRegression(swlnoise, swrnoise, REJnoise, LEJnoise, RSJnoise, ...
    LSJnoise, x_wrnoise, x_wlnoise, dx_wrnoise, dx_wlnoise)

% saveas(gcf,['plotresults/data_trial_extra_noise',num2str(idfig),'.pdf'])

set(gcf,'renderer','Painters')
% saveas(gcf,['plotresults/gpr_trial',num2str(idfig),'.pdf'])
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,['plotresults/data_trial_extra_noise',num2str(idfig)],'-dpdf','-r0')

function out = addNoise(sig, amplitude)

    if isa(sig,'struct')
        fnames = fieldnames(sig);
        for i=1:8
            sig.(fnames{i}) = sig.(fnames{i}) + amplitude*rand(size(sig.(fnames{i})));
        end
        out = sig;
    else
        out = sig + amplitude*rand(size(sig));
    end
end

