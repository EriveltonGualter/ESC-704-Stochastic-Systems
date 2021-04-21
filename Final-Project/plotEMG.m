function plotEMG(emg)

    fig = figure; fig.Renderer = 'painters';
    for i=2:11
        subplot(5,2,i-1); plot(emg(:,1), emg(:,i)); 
        ylabel(['EMG', num2str(i-1)]);
    end
end