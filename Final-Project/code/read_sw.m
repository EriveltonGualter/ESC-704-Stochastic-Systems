function [ sw ] = read_sw(raw)   
% sw = read_sw(raw)   
% 
%     This function creates the structure data for SmartWheel data
% 
%     INPUTS:
%         raw = [n, 24] = Raw Data 
%
%     OUTPUTS:
%         sw = Forces and Moments structure
%             .Sn = Sample Number 
%             .av = 1/20 second running average velocity [m/s]
%             .Fx = Force [N]
%             .Fy = Force [N]
%             .Fz = Force [N]
%             .Mx = Moment [N.m]
%             .My = Moment [N.m]
%             .Mz = Moment [N.m]

    for i=1:length(raw)
        sw.Sn(i,1) = raw(i,1);
        sw.av(i,1) = raw(i,5);
        sw.Fx(i,1) = raw(i,19);
        sw.Fy(i,1) = raw(i,20);
        sw.Fz(i,1) = raw(i,21);
        sw.Mx(i,1) = raw(i,22);
        sw.My(i,1) = raw(i,23);
        sw.Mz(i,1) = raw(i,24);
    end
end