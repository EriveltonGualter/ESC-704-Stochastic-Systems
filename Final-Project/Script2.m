clear all

run('Script1_Preprocessing.m')
close all

figure; hold on;
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

view([45 30]); axis equal; axis([-1.5e3 5e3 -.5e3 .5e3 0 1.5e3]); xlabel('x'); ylabel('y'); zlabel('z');

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
PM(any(isnan(PM), 2), :) = [];

t = PM(:,end)-PM(1,end);

j = 1;
for i=1:28
    M(i).x = PM(:,j); j=j+1;
    M(i).y = PM(:,j); j=j+1;
    M(i).z = PM(:,j); j=j+1;
end

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

cla; hold on; box on
clearvars -global;
clear global;

view([62 17]);
% view(2);
global h;

varnames = propulsion.Properties.VariableNames;

cp = [0 0 0];

i = 1;
time = 0;
tic 
while time < t(end)
% for i=1:size(PM,1)
    
    if isempty(h)
        [px,py,pz] = getSegment(MK.IJ,  MK.PX,t,time);  h(1) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.IJ,  MK.RAC,t,time); h(2) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.RAC, MK.RLE,t,time); h(3) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.RLE, MK.RRS,t,time); h(4) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.IJ,  MK.RAC,t,time); h(5) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.RAC, MK.RLE,t,time); h(6) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        [px,py,pz] = getSegment(MK.RLE, MK.RRS,t,time); h(7) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
        
        [px,py,pz] = getSegment(MK.RW, MK.LW,t,time);   h(8) = plot3(px,py,pz,'LineWidth',1,'Color',cp);
    else
        [px,py,pz] = getSegment(MK.IJ,  MK.PX,t,time);  set(h(1),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.IJ,  MK.RAC,t,time); set(h(2),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.RAC, MK.RLE,t,time); set(h(3),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.RLE, MK.RRS,t,time); set(h(4),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.IJ,  MK.LAC,t,time); set(h(5),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.LAC, MK.LLE,t,time); set(h(6),'xData',px,'yData',py,'zData',pz);
        [px,py,pz] = getSegment(MK.LLE, MK.LRS,t,time); set(h(7),'xData',px,'yData',py,'zData',pz);
        
        [px,py,pz] = getSegment(MK.RW, MK.LW,t,time);   set(h(8),'xData',px,'yData',py,'zData',pz);
        
    end
    
%     view([45 30]); 
    axis equal; axis([-1.5e3 5e3 -.5e3 .5e3 0 1.5e3]); 
    set(gca,'xtick',[],'ytick',[],'ztick',[])
%     xlabel('x'); ylabel('y'); zlabel('z'); 

    drawnow
    
    % Update current time
    time = toc;
end


function [px,py,pz] = getSegment(A,B,t,time)   
    max = interp1(t',A.x',time')';
    may = interp1(t',A.y',time')';
    maz = interp1(t',A.z',time')';
    
    mbx = interp1(t',B.x',time')';
    mby = interp1(t',B.y',time')';
    mbz = interp1(t',B.z',time')';
    
    px = [mbx max];
    py = [mby may];
    pz = [mbz maz];
end

