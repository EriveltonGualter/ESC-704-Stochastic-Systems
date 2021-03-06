function [ emg ] = read_emg(raw)
    % emg = read_emg(raw)   
% 
%     This function creates the structure data for EMG data
% 
%     INPUTS:
%         raw = [n, 24] = Raw Data 
%
%     OUTPUTS:
%         sw = Forces and Moments structure
%             .Sn = Sample Number 
%             .av = 1/20 second running average velocity [m/s]
%             .Fx = Force [N]

    % SSynchronization of arrays raw.data and raw.textdat
    s = size(raw.data);
    raw.data = cat(1, ones(1,s(1,2)), raw.data);
    
    ii = 1;
    for i=1:length(raw.data)
        if raw.textdata(i,1) == "EMG"   % Just save the EMG data
            
            emg.frame(ii,1) = raw.data(i,1);
            emg.t(ii,1) = raw.data(i,2);
            emg.ch1(ii,1) = raw.data(i,3);
            emg.ch2(ii,1) = raw.data(i,4);
            emg.ch3(ii,1) = raw.data(i,5);
            emg.ch4(ii,1) = raw.data(i,6);
            emg.ch5(ii,1) = raw.data(i,7);
            emg.ch6(ii,1) = raw.data(i,8);
            emg.ch7(ii,1) = raw.data(i,9);
            emg.ch8(ii,1) = raw.data(ii,10);
            emg.ch9(ii,1) = raw.data(ii,11);
            emg.ch10(ii,1) = raw.data(ii,12);   
            
            ii = ii + 1;
        end
    end
end