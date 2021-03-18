function plots(swl, swr )
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
       
    subplot(421); hold on; box on; ylabel('Fx [N])'); title('Forces applied in the left wheel')
    plot(tswl, swl.Fx, 'LineWidth', 1, 'Color', 'b');
    plot(tswr, swr.Fx, 'LineWidth', 1, 'Color', 'r');  
                   
    subplot(423); hold on; box on; ylabel('Fy [N])');
    plot(tswl, swl.Fy, 'LineWidth', 1, 'Color', 'b');
    plot(tswr, swr.Fy, 'LineWidth', 1, 'Color', 'r');  
    
    subplot(425); hold on; box on; ylabel('Fz [N])');
    plot(tswl, swl.Fz, 'LineWidth', 1, 'Color', 'b');
    plot(tswr, swr.Fz, 'LineWidth', 1, 'Color', 'r');  
    
    subplot(427); hold on; box on; ylabel('Fr [N])');
    plot(tswl, sqrt(swl.Fx.*swl.Fx + swl.Fy.*swl.Fy + swl.Fz.*swl.Fz), 'LineWidth', 1, 'Color', 'b');
    plot(tswr, sqrt(swr.Fx.*swr.Fx + swr.Fy.*swr.Fy + swr.Fz.*swr.Fz), 'LineWidth', 1, 'Color', 'r');    
    xlabel(strXlabel) 
      
    
    %% Plot Moments-left x Sample
       
    subplot(422); hold on; box on; ylabel('Mx [Nm])'); title('Moments applied in the left wheel')
    plot(tswl, swl.Mx, 'LineWidth', 1, 'Color', 'b');
    plot(tswr, swr.Mx, 'LineWidth', 1, 'Color', 'r');      

    subplot(424); hold on; box on; ylabel('My [Nm])');
    plot(tswl, swl.My, 'LineWidth', 1, 'Color', 'b');
    plot(tswr, swr.My, 'LineWidth', 1, 'Color', 'r');    
    
    subplot(426); hold on; box on; ylabel('Mz [Nm])');
    plot(tswl, swl.Mz, 'LineWidth', 1, 'Color', 'b');
    plot(tswr, swr.Mz, 'LineWidth', 1, 'Color', 'r');      
    
    subplot(428); hold on; box on; ylabel('Mr [Nm])');
    plot(tswl, sqrt(swl.Mx.*swl.Mx + swl.My.*swl.My + swl.Mz.*swl.Mz), 'LineWidth', 1, 'Color', 'b');
    plot(tswr, sqrt(swr.Mx.*swr.Mx + swr.My.*swr.My + swr.Mz.*swr.Mz), 'LineWidth', 1, 'Color', 'r');    
    xlabel(strXlabel) 
end

