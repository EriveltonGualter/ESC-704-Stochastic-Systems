% clear all

addpath('helpers');

% trial = 1;
run('Script1_Preprocessing.m');

%%
set(0, 'CurrentFigure', fig(1))

rgbl = [0,0,0,.2];
plot3(propulsion.IJ,	propulsion.VarName5	,	propulsion.VarName4,    'Color',rgbl);
plot3(propulsion.PX,	propulsion.VarName8	,	propulsion.VarName7,    'Color',rgbl);
plot3(propulsion.C7,	propulsion.VarName11,	propulsion.VarName10,   'Color',rgbl);
plot3(propulsion.T8,	propulsion.VarName14,	propulsion.VarName13,   'Color',rgbl);
plot3(propulsion.RAC,	propulsion.VarName17,	propulsion.VarName16,   'Color',rgbl);
plot3(propulsion.RAS,	propulsion.VarName20,	propulsion.VarName19,   'Color',rgbl);
plot3(propulsion.RPS,	propulsion.VarName23,	propulsion.VarName22,   'Color',rgbl);
plot3(propulsion.RAI,	propulsion.VarName26,	propulsion.VarName25,   'Color',rgbl);
plot3(propulsion.RPI,	propulsion.VarName29,	propulsion.VarName28,   'Color',rgbl);
plot3(propulsion.RLE,	propulsion.VarName32,	propulsion.VarName31,   'Color',rgbl);
plot3(propulsion.RME,	propulsion.VarName35,	propulsion.VarName34,   'Color',rgbl);
plot3(propulsion.RRS,	propulsion.VarName38,	propulsion.VarName37,   'Color',rgbl);
plot3(propulsion.RUS,	propulsion.VarName41,	propulsion.VarName40,   'Color',rgbl);
plot3(propulsion.R3M,	propulsion.VarName44,	propulsion.VarName43,   'Color',rgbl);
plot3(propulsion.R5M,	propulsion.VarName47,	propulsion.VarName46,   'Color',rgbl);
plot3(propulsion.LAC,	propulsion.VarName50,	propulsion.VarName49,   'Color',rgbl);
plot3(propulsion.LAS,	propulsion.VarName53,	propulsion.VarName52,   'Color',rgbl);
plot3(propulsion.LPS,	propulsion.VarName56,	propulsion.VarName55,   'Color',rgbl);
plot3(propulsion.LAI,	propulsion.VarName59,	propulsion.VarName58,   'Color',rgbl);
plot3(propulsion.LPI,	propulsion.VarName62,	propulsion.VarName61,   'Color',rgbl);
plot3(propulsion.LLE,	propulsion.VarName65,	propulsion.VarName64,   'Color',rgbl);
plot3(propulsion.LME,	propulsion.VarName68,	propulsion.VarName67,   'Color',rgbl);
plot3(propulsion.LRS,	propulsion.VarName71,	propulsion.VarName70,   'Color',rgbl);
plot3(propulsion.LUS,	propulsion.VarName74,	propulsion.VarName73,   'Color',rgbl);
plot3(propulsion.L3M,	propulsion.VarName77,	propulsion.VarName76,   'Color',rgbl);
plot3(propulsion.L5M,	propulsion.VarName80,	propulsion.VarName79,   'Color',rgbl);
plot3(propulsion.RW,	propulsion.VarName83,	propulsion.VarName82,   'Color',rgbl);
plot3(propulsion.LW,	propulsion.VarName86,	propulsion.VarName85,   'Color',rgbl);

view([62 17]);
% view([45 30]);
axis equal; axis([-1.5e3 5e3 -.5e3 .5e3 0 1.5e3]); 
% xlabel('x'); ylabel('y'); zlabel('z');

%%
PosMarker = [ ...
    propulsion.IJ,	propulsion.VarName5	,	propulsion.VarName4, ...
    propulsion.PX,	propulsion.VarName8	,	propulsion.VarName7, ...
    propulsion.C7,	propulsion.VarName11,	propulsion.VarName10, ...
    propulsion.T8,	propulsion.VarName14,	propulsion.VarName13, ...
    propulsion.RAC,	propulsion.VarName17,	propulsion.VarName16, ...
    propulsion.RAS,	propulsion.VarName20,	propulsion.VarName19, ...
    propulsion.RPS,	propulsion.VarName23,	propulsion.VarName22, ...
    propulsion.RAI,	propulsion.VarName26,	propulsion.VarName25, ...
    propulsion.RPI,	propulsion.VarName29,	propulsion.VarName28, ...
    propulsion.RLE,	propulsion.VarName32,	propulsion.VarName31, ...
    propulsion.RME,	propulsion.VarName35,	propulsion.VarName34, ...
    propulsion.RRS,	propulsion.VarName38,	propulsion.VarName37, ...
    propulsion.RUS,	propulsion.VarName41,	propulsion.VarName40, ...
    propulsion.R3M,	propulsion.VarName44,	propulsion.VarName43, ...
    propulsion.R5M,	propulsion.VarName47,	propulsion.VarName46, ...
    propulsion.LAC,	propulsion.VarName50,	propulsion.VarName49, ...
    propulsion.LAS,	propulsion.VarName53,	propulsion.VarName52, ...
    propulsion.LPS,	propulsion.VarName56,	propulsion.VarName55, ...
    propulsion.LAI,	propulsion.VarName59,	propulsion.VarName58, ...
    propulsion.LPI,	propulsion.VarName62,	propulsion.VarName61, ...
    propulsion.LLE,	propulsion.VarName65,	propulsion.VarName64, ...
    propulsion.LME,	propulsion.VarName68,	propulsion.VarName67, ...
    propulsion.LRS,	propulsion.VarName71,	propulsion.VarName70, ...
    propulsion.LUS,	propulsion.VarName74,	propulsion.VarName73, ...
    propulsion.L3M,	propulsion.VarName77,	propulsion.VarName76, ...
    propulsion.L5M,	propulsion.VarName80,	propulsion.VarName79, ...
    propulsion.RW,	propulsion.VarName83,	propulsion.VarName82, ...
    propulsion.LW,	propulsion.VarName86,	propulsion.VarName85, ...
    propulsion.Time];

PM = PosMarker;
PM(1:find(PM(:,end)==0, 1, 'last'), :) = [];

PM(1:find(PM(:,end)<=tanc(1), 1, 'last'), :) = [];
PM(find(PM(:,end)<=tanc(end), 1, 'last'):end, :) = [];

t = PM(:,end)-PM(1,end);

% PM(any(isnan(PM), 2), :) = [];

%%
f1 = length(t)/t(end);
f2 = length(swl.Fx)/t(end);

t2 = 0:1/f2:(length(swl.Fx)-1)*1/f2;
t2new = 0:1/f1:t2(end);

% PM(find(PM(:,end)<=t2new(end), 1, 'last'):end, :) = [];

swl2.Fx = interp1(t2,swl.Fx,t2new);
swl2.Fy = interp1(t2,swl.Fy,t2new);
swl2.Fz = interp1(t2,swl.Fz,t2new);
swl2.F  = sqrt(swl2.Fx.*swl2.Fx + swl2.Fy.*swl2.Fy + swl2.Fz.*swl2.Fz);
swl2.Mx = interp1(t2,swl.Mx,t2new);
swl2.My = interp1(t2,swl.My,t2new);
swl2.Mz = interp1(t2,swl.Mz,t2new);
swl2.M  = sqrt(swl2.Mx.*swl2.Mx + swl2.My.*swl2.My + swl2.Mz.*swl2.Mz);
swl2.t  = t2new;

t2 = 0:1/f2:(length(swr.Fx)-1)*1/f2;
t2new = 0:1/f1:t2(end);

swr2.Fx = interp1(t2,swr.Fx,t2new);
swr2.Fy = interp1(t2,swr.Fy,t2new);
swr2.Fz = interp1(t2,swr.Fz,t2new);
swr2.F  = sqrt(swr2.Fx.*swr2.Fx + swr2.Fy.*swr2.Fy + swr2.Fz.*swr2.Fz);
swr2.Mx = interp1(t2,swr.Mx,t2new);
swr2.My = interp1(t2,swr.My,t2new);
swr2.Mz = interp1(t2,swr.Mz,t2new);
swr2.M  = sqrt(swr2.Mx.*swr2.Mx + swr2.My.*swr2.My + swr2.Mz.*swr2.Mz);
swr2.t  = t2new;

set(0, 'CurrentFigure', fig(2))
plots( swl2, swr2 );

%%
idx = find((any(isnan(PM), 2)) ~= 1);

swl2.Fx = swl2.Fx(idx);
swl2.Fy = swl2.Fy(idx);
swl2.Fz = swl2.Fz(idx);
swl2.F  = swl2.F(idx);
swl2.Mx = swl2.Mx(idx);
swl2.My = swl2.My(idx);
swl2.Mz = swl2.Mz(idx);
swl2.M  = swl2.M(idx);
swl2.t  = swl2.t(idx);

swr2.Fx = swr2.Fx(idx);
swr2.Fy = swr2.Fy(idx);
swr2.Fz = swr2.Fz(idx);
swr2.F  = swr2.F(idx);
swr2.Mx = swr2.Mx(idx);
swr2.My = swr2.My(idx);
swr2.Mz = swr2.Mz(idx);
swr2.M  = swr2.M(idx);
swr2.t  = swr2.t(idx);

PM = PM(idx,:);
t = t(idx);
t = t - t(1);

%%
j = 1;
for i=1:28
    M(i).x = PM(:,j); j=j+1;
    M(i).y = PM(:,j); j=j+1;
    M(i).z = PM(:,j); j=j+1;
end

%%

MK.IJ.x  = M(1).x;  MK.IJ.y  = M(1).y;  MK.IJ.z  = M(1).z; 
MK.PX.x  = M(2).x;  MK.PX.y  = M(2).y;  MK.PX.z  = M(2).z; 
MK.C7.x  = M(3).x;  MK.C7.y  = M(3).y;  MK.C7.z  = M(3).z; 
MK.T8.x  = M(4).x;  MK.T8.y  = M(4).y;  MK.T8.z  = M(4).z; 
MK.RAC.x = M(5).x;  MK.RAC.y = M(5).y;  MK.RAC.z = M(5).z; 
MK.RAS.x = M(6).x;  MK.RAS.y = M(6).y;  MK.RAS.z = M(6).z; 
MK.RPS.x = M(7).x;  MK.RPS.y = M(7).y;  MK.RPS.z = M(7).z; 
MK.RAI.x = M(8).x;  MK.RAI.y = M(8).y;  MK.RAI.z = M(8).z; 
MK.RPI.x = M(9).x;  MK.RPI.y = M(9).y;  MK.RPI.z = M(9).z; 
MK.RLE.x = M(10).x; MK.RLE.y = M(10).y; MK.RLE.z = M(10).z; 
MK.RME.x = M(11).x; MK.RME.y = M(11).y; MK.RME.z = M(11).z; 
MK.RRS.x = M(12).x; MK.RRS.y = M(12).y; MK.RRS.z = M(12).z; 
MK.RUS.x = M(13).x; MK.RUS.y = M(13).y; MK.RUS.z = M(13).z; 
MK.R3M.x = M(14).x; MK.R3M.y = M(14).y; MK.R3M.z = M(14).z; 
MK.R5M.x = M(15).x; MK.R5M.y = M(15).y; MK.R5M.z = M(15).z; 
MK.LAC.x = M(16).x; MK.LAC.y = M(16).y; MK.LAC.z = M(16).z; 
MK.LAS.x = M(17).x; MK.LAS.y = M(17).y; MK.LAS.z = M(17).z; 
MK.LPS.x = M(18).x; MK.LPS.y = M(18).y; MK.LPS.z = M(18).z; 
MK.LAI.x = M(19).x; MK.LAI.y = M(19).y; MK.LAI.z = M(19).z; 
MK.LPI.x = M(20).x; MK.LPI.y = M(20).y; MK.LPI.z = M(20).z; 
MK.LLE.x = M(21).x; MK.LLE.y = M(21).y; MK.LLE.z = M(21).z; 
MK.LME.x = M(22).x; MK.LME.y = M(22).y; MK.LME.z = M(22).z; 
MK.LRS.x = M(23).x; MK.LRS.y = M(23).y; MK.LRS.z = M(23).z; 
MK.LUS.x = M(24).x; MK.LUS.y = M(24).y; MK.LUS.z = M(24).z; 
MK.L3M.x = M(25).x; MK.L3M.y = M(25).y; MK.L3M.z = M(25).z; 
MK.L5M.x = M(26).x; MK.L5M.y = M(26).y; MK.L5M.z = M(26).z; 
MK.RW.x  = M(27).x; MK.RW.y  = M(27).y; MK.RW.z  = M(27).z; 
MK.LW.x  = M(28).x; MK.LW.y  = M(28).y; MK.LW.z  = M(28).z; 

%%
set(0, 'CurrentFigure', fig(3))
cla; hold on; box on
clearvars -global;
clear global;

% view([62 17]);
% view(2);
view([0 0]);
view([0 0 1])
global h hw;

varnames = propulsion.Properties.VariableNames;

cp = [0 0 0];

simRealTime = 3;
switch simRealTime
    case 1
        N = t(end);
        i = 0;
        time = 0;
    case 2
        N = length(MK.C7.x);
        i = 1;
    case 3
        N = length(MK.C7.x);
        time = Inf;
end

tic 
while time < N
    
    if isempty(h)
%         [px,py,pz] = getSegment(MK.T8,  MK.C7,t,time,i);  h(1) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.IJ,  MK.PX,t,time,i);  h(2) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.IJ,  MK.RAC,t,time,i); h(3) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.RAC, MK.RLE,t,time,i); h(4) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.RLE, MK.RRS,t,time,i); h(5) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.IJ,  MK.LAC,t,time,i); h(6) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.LAC, MK.LLE,t,time,i); h(7) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.LLE, MK.LRS,t,time,i); h(8) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        
        [px,py,pz] = getSegment(MK.RW, MK.LW,t,time,i);   h(9) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
    else
%         [px,py,pz] = getSegment(MK.T8,  MK.C7,t,time,i);  set(h(1),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.IJ,  MK.PX,t,time,i);  set(h(2),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.IJ,  MK.RAC,t,time,i); set(h(3),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.RAC, MK.RLE,t,time,i); set(h(4),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.RLE, MK.RRS,t,time,i); set(h(5),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.IJ,  MK.LAC,t,time,i); set(h(6),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.LAC, MK.LLE,t,time,i); set(h(7),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.LLE, MK.LRS,t,time,i); set(h(8),'xData',px,'yData',py,'zData',pz);
        
        [px,py,pz] = getSegment(MK.RW, MK.LW,t,time,i);   set(h(9),'xData',px,'yData',py,'zData',pz);
        
    end
    
    hw = plotwheels(MK.RW,MK.LW,t,time,i,hw);
    
    view([45 30]); 
    axis equal; axis([-1.5e3 5e3 -.5e3 .5e3 0 1.5e3]); 
    set(gca,'xtick',[],'ytick',[],'ztick',[])
%     xlabel('x'); ylabel('y'); zlabel('z'); 

    drawnow;
    
    % Update current time
    if simRealTime == 1
        time = toc;
    else
        time = i + 1;
        i = time;
    end
end
%%
for i=1:length(MK.RAC.x)
    P1 = MK.RAC;
    P2 = MK.RLE;
    P3 = MK.RLE;
    P4 = MK.RRS;
    REJ(i) = getAngle(P1, P2, P3, P4, i);   % Right Elbow Joint 
    
    P1 = MK.LAC;
    P2 = MK.LLE;
    P3 = MK.LLE;
    P4 = MK.LRS;
    LEJ(i) = getAngle(P1, P2, P3, P4, i);   % Left Elbow Joint 

    P1 = MK.RAC;
    P2 = MK.RLE;
    P3 = MK.PX;
    P4 = MK.IJ;
    RSJ(i) = getAngle(P1, P2, P3, P4, i);   % Right Shoulder Joint 
    
    P1 = MK.LAC;
    P2 = MK.LLE;
    P3 = MK.PX;
    P4 = MK.IJ;
    LSJ(i) = getAngle(P1, P2, P3, P4, i);   % LEFT Shoulder Joint 
end

set(0, 'CurrentFigure', fig(4))
subplot(221); box on; hold on; plot(t, rad2deg(REJ));
subplot(222); box on; hold on; plot(t, rad2deg(LEJ));
subplot(223); box on; hold on; plot(t, rad2deg(RSJ));
subplot(224); box on; hold on; plot(t, rad2deg(LSJ));

%%

set(0, 'CurrentFigure', fig(5))

subplot(221); hold on; box on; peaks = plotStride(REJ);
subplot(222); hold on; box on; plotStride(LEJ,peaks);
subplot(223); hold on; box on; plotStride(RSJ,peaks);
subplot(224); hold on; box on; plotStride(LSJ,peaks);

%%
for i=1:length(fig)
    fig(i).Visible = 'off';
end



