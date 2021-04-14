function plotDTforRegression(swl, swr, REJ2, LEJ2, RSJ, LSJ, x_wr, x_wl, dx_wr, dx_wl)
% sw = read_sw(raw)   
% 
%     This function is ...
% 
%     INPUTS:
%         trc = [1, 1] = Makers Position
%         swl = [1, 1] = SmarWheel Data - left
%         swr = [1, 1] = SmarWheel Data - right

    if isfield(swl,'t')
        tswl = swl.t;
        tswr = swr.t;
        strXlabel = 'Time, s';
    else
        tswl = swl.Sn;
        tswr = swr.Sn;
        strXlabel = 'Sample';
    end

    %% Plot Forces-left x Sample

    subplot(521); hold on; box on; plot(tswr, x_wl); ylabel('m'); title('Displacement Left'); 
    subplot(522); hold on; box on; plot(tswl, x_wr); ylabel('m'); title('Displacement Right'); 
    subplot(523); hold on; box on; plot(tswr, -REJ2(:,1)); ylabel('Angle'); title('Alpha Left'); 
    subplot(524); hold on; box on; plot(tswl, -LEJ2(:,1)); ylabel('Angle'); title('Alpha Right'); 
    subplot(525); hold on; box on; plot(tswr, 1.8*RSJ); ylabel('Angle'); title('Beta Left')
    subplot(526); hold on; box on; plot(tswl, 1.8*LSJ); ylabel('Angle'); title('Beta Right')

    subplot(527); hold on; box on; ylabel('Force [N]'); title('Forces applied in the left wheel')
    plot(tswl, swl.Fx, 'LineWidth', 1);
    plot(tswl, swl.Fy, 'LineWidth', 1);
    plot(tswl, swl.Fz, 'LineWidth', 1);
    plot(tswl, sqrt(swl.Fx.*swl.Fx + swl.Fy.*swl.Fy + swl.Fz.*swl.Fz), 'LineWidth', 1);
    
    subplot(528); hold on; box on; ylabel('Force [N]'); title('Forces applied in the right wheel')
    plot(tswr, swr.Fx, 'LineWidth', 1);  
    plot(tswr, swr.Fy, 'LineWidth', 1);  
    plot(tswr, swr.Fz, 'LineWidth', 1);  
    plot(tswr, sqrt(swr.Fx.*swr.Fx + swr.Fy.*swr.Fy + swr.Fz.*swr.Fz), 'LineWidth', 1);    
    xlabel(strXlabel) 
    
    subplot(529); hold on; box on; ylabel('Moment [Nm]'); title('Moments applied in the left wheel')
    plot(tswl, swl.Mx, 'LineWidth', 1);
    plot(tswl, swl.My, 'LineWidth', 1);
    plot(tswl, swl.Mz, 'LineWidth', 1);
    plot(tswl, sqrt(swl.Mx.*swl.Mx + swl.My.*swl.My + swl.Mz.*swl.Mz), 'LineWidth', 1);
    
    subplot(5,2,10); hold on; box on; ylabel('Moment [Nm]'); title('Moments applied in the right wheel')
    plot(tswr, swr.Mx, 'LineWidth', 1);      
    plot(tswr, swr.My, 'LineWidth', 1);    
    plot(tswr, swr.Mz, 'LineWidth', 1);      
    plot(tswr, sqrt(swr.Mx.*swr.Mx + swr.My.*swr.My + swr.Mz.*swr.Mz), 'LineWidth', 1);    
    
    
end

