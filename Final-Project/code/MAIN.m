clear, clc

addpath('data')

% Load Raw Data
trc_raw = dlmread('wc2.trc', '', 6, 0); % Raw data from Tracker Position
emg_raw = importdata('wc2-Delsys.csv', ',', 1,1); % Raw data from EMG
swl_raw = dlmread('wc_SW1.csv', ',');   % Raw data from Left SmartWheel
swr_raw = dlmread('wc_SW2.csv', ',');   % Raw data from Right SmartWheel

% Cleaning and Structuring Data
trc = read_trc( trc_raw );  % structure tracker data
emg = read_emg( emg_raw );  % structure EMG data
swl = read_sw( swl_raw );   % strucutre of left wheel 
swr = read_sw( swr_raw );   % strucutre of right wheel 
