clear all
close all
clc

for i=1:5
    fig(i) = figure; fig(i).Renderer = 'painters'; hold on; box on; fig(i).Visible = 'off';
end

for trial=2:4
    clearvars all -except trial fig
    disp(trial)
    Script2
end

for i=1:length(fig)
    fig(i).Visible = 'on';
end
