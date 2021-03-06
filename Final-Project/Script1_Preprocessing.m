addpath(genpath(pwd))

if ~exist('trial','var')
    trial = 2;
end


%% Import data from text file
% Script for importing data from the following text file:


%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 102);

% Specify range and delimiter
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Var1", "Time", "IJ", "VarName4", "VarName5", "PX", "VarName7", "VarName8", "C7", "VarName10", "VarName11", "T8", "VarName13", "VarName14", "RAC", "VarName16", "VarName17", "RAS", "VarName19", "VarName20", "RPS", "VarName22", "VarName23", "RAI", "VarName25", "VarName26", "RPI", "VarName28", "VarName29", "RLE", "VarName31", "VarName32", "RME", "VarName34", "VarName35", "RRS", "VarName37", "VarName38", "RUS", "VarName40", "VarName41", "R3M", "VarName43", "VarName44", "R5M", "VarName46", "VarName47", "LAC", "VarName49", "VarName50", "LAS", "VarName52", "VarName53", "LPS", "VarName55", "VarName56", "LAI", "VarName58", "VarName59", "LPI", "VarName61", "VarName62", "LLE", "VarName64", "VarName65", "LME", "VarName67", "VarName68", "LRS", "VarName70", "VarName71", "LUS", "VarName73", "VarName74", "L3M", "VarName76", "VarName77", "L5M", "VarName79", "VarName80", "RW", "VarName82", "VarName83", "LW", "VarName85", "VarName86", "Var87", "Var88", "Var89", "Var90", "Var91", "Var92", "Var93", "Var94", "Var95", "Var96", "Var97", "Var98", "Var99", "Var100", "Var101", "Var102"];
opts.SelectedVariableNames = ["Time", "IJ", "VarName4", "VarName5", "PX", "VarName7", "VarName8", "C7", "VarName10", "VarName11", "T8", "VarName13", "VarName14", "RAC", "VarName16", "VarName17", "RAS", "VarName19", "VarName20", "RPS", "VarName22", "VarName23", "RAI", "VarName25", "VarName26", "RPI", "VarName28", "VarName29", "RLE", "VarName31", "VarName32", "RME", "VarName34", "VarName35", "RRS", "VarName37", "VarName38", "RUS", "VarName40", "VarName41", "R3M", "VarName43", "VarName44", "R5M", "VarName46", "VarName47", "LAC", "VarName49", "VarName50", "LAS", "VarName52", "VarName53", "LPS", "VarName55", "VarName56", "LAI", "VarName58", "VarName59", "LPI", "VarName61", "VarName62", "LLE", "VarName64", "VarName65", "LME", "VarName67", "VarName68", "LRS", "VarName70", "VarName71", "LUS", "VarName73", "VarName74", "L3M", "VarName76", "VarName77", "L5M", "VarName79", "VarName80", "RW", "VarName82", "VarName83", "LW", "VarName85", "VarName86"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var87", "Var88", "Var89", "Var90", "Var91", "Var92", "Var93", "Var94", "Var95", "Var96", "Var97", "Var98", "Var99", "Var100", "Var101", "Var102"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var87", "Var88", "Var89", "Var90", "Var91", "Var92", "Var93", "Var94", "Var95", "Var96", "Var97", "Var98", "Var99", "Var100", "Var101", "Var102"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, ["L5M", "VarName79", "VarName80"], "TrimNonNumeric", true);
opts = setvaropts(opts, ["L5M", "VarName79", "VarName80"], "ThousandsSeparator", ",");

% Import the data
if trial < 10
    propulsion = readtable(append("trial0",num2str(trial),"/wc",num2str(trial),".trc"), opts);
else
    propulsion = readtable(append("trial",num2str(trial),"/wc",num2str(trial),".trc"), opts);
end

%% Clear temporary variables
clear opts

%%
% Load Raw Data
trc_raw = dlmread(      append('wc',    num2str(trial),'.trc'), '', 6, 0);          % Raw data from Tracker Position
emg_raw = importdata(   append('wc',    num2str(trial),'-Delsys 1.csv'), ',', 1,1); % Raw data from EMG
if trial < 10 
    swl_raw = dlmread(      append('SW1-0', num2str(trial),'.csv'), ',');           % Raw data from Left SmartWheel
    swr_raw = dlmread(      append('SW2-0', num2str(trial),'.csv'), ',');           % Raw data from Right SmartWheel
else
    swl_raw = dlmread(      append('SW1-', num2str(trial),'.csv'), ',');           % Raw data from Left SmartWheel
    swr_raw = dlmread(      append('SW2-', num2str(trial),'.csv'), ',');           % Raw data from Right SmartWheel
end

% Cleaning and Structuring Data
trc = read_trc( trc_raw );  % structure tracker data
swl = read_sw( swl_raw );   % strucutre of left wheel 
swr = read_sw( swr_raw );   % strucutre of right wheel 
emg = read_emg( trial );    % structure EMG data
tanc= read_anc( trial );

%%
% figure; hold on;
% plot3(propulsion.IJ,	propulsion.VarName5	,	propulsion.VarName4);
% plot3(propulsion.PX,	propulsion.VarName8	,	propulsion.VarName7);
% plot3(propulsion.C7,	propulsion.VarName11,	propulsion.VarName10);
% plot3(propulsion.T8,	propulsion.VarName14,	propulsion.VarName13);
% plot3(propulsion.RAC,	propulsion.VarName17,	propulsion.VarName16);
% plot3(propulsion.RAS,	propulsion.VarName20,	propulsion.VarName19);
% plot3(propulsion.RPS,	propulsion.VarName23,	propulsion.VarName22);
% plot3(propulsion.RAI,	propulsion.VarName26,	propulsion.VarName25);
% plot3(propulsion.RPI,	propulsion.VarName29,	propulsion.VarName28);
% plot3(propulsion.RLE,	propulsion.VarName32,	propulsion.VarName31);
% plot3(propulsion.RME,	propulsion.VarName35,	propulsion.VarName34);
% plot3(propulsion.RRS,	propulsion.VarName38,	propulsion.VarName37);
% plot3(propulsion.RUS,	propulsion.VarName41,	propulsion.VarName40);
% plot3(propulsion.R3M,	propulsion.VarName44,	propulsion.VarName43);
% plot3(propulsion.R5M,	propulsion.VarName47,	propulsion.VarName46);
% plot3(propulsion.LAC,	propulsion.VarName50,	propulsion.VarName49);
% plot3(propulsion.LAS,	propulsion.VarName53,	propulsion.VarName52);
% plot3(propulsion.LPS,	propulsion.VarName56,	propulsion.VarName55);
% plot3(propulsion.LAI,	propulsion.VarName59,	propulsion.VarName58);
% plot3(propulsion.LPI,	propulsion.VarName62,	propulsion.VarName61);
% plot3(propulsion.LLE,	propulsion.VarName65,	propulsion.VarName64);
% plot3(propulsion.LME,	propulsion.VarName68,	propulsion.VarName67);
% plot3(propulsion.LRS,	propulsion.VarName71,	propulsion.VarName70);
% plot3(propulsion.LUS,	propulsion.VarName74,	propulsion.VarName73);
% plot3(propulsion.L3M,	propulsion.VarName77,	propulsion.VarName76);
% plot3(propulsion.L5M,	propulsion.VarName80,	propulsion.VarName79);
% plot3(propulsion.RW,	propulsion.VarName83,	propulsion.VarName82);
% plot3(propulsion.LW,	propulsion.VarName86,	propulsion.VarName85);
% 
% view([45 30]); axis equal; axis([-1.5e3 5e3 -.5e3 .5e3 0 1.5e3]); xlabel('x'); ylabel('y'); zlabel('z');

