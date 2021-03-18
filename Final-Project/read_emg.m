function [ emg ] = read_emg(trial)

    %% Set up the Import Options and import the data
    opts = delimitedTextImportOptions("NumVariables", 20);

    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["Var1", "Signal", "Frame", "Seconds", "Channels_1", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20"];
    opts.SelectedVariableNames = ["Signal", "Frame", "Seconds", "Channels_1", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14"];
    opts.VariableTypes = ["string", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string"];

    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Specify variable properties
    opts = setvaropts(opts, ["Var1", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20"], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, ["Var1", "Signal", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20"], "EmptyFieldRule", "auto");

    % Import the data
    if trial < 10
        emg_raw = readtable(append("trial0",num2str(trial),"/wc",num2str(trial),"-Delsys 1.csv"), opts);
    else
        emg_raw = readtable(append("trial",num2str(trial),"/wc",num2str(trial),"-Delsys 1.csv"), opts);
    end

    idx = find(emg_raw.Signal == 'EMG'); 
    
    emg(:,1)    = emg_raw.Seconds(idx);
    emg(:,2)    = emg_raw.Channels_1(idx);
    emg(:,3)    = emg_raw.VarName6(idx);
    emg(:,4)    = emg_raw.VarName7(idx);
    emg(:,5)    = emg_raw.VarName8(idx);
    emg(:,6)    = emg_raw.VarName9(idx);
    emg(:,7)    = emg_raw.VarName10(idx);
    emg(:,8)    = emg_raw.VarName11(idx);
    emg(:,9)    = emg_raw.VarName12(idx);
    emg(:,10)   = emg_raw.VarName13(idx);
    emg(:,11)   = emg_raw.VarName14(idx);