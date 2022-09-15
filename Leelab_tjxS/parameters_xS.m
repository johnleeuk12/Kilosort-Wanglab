function x = parameters_xS(PC, username,fpath2)

%% Do not change
x.PC_name = PC;
x.animal_name = username;
[x.fpath, x.savepath] = directories(x.PC_name);
x.fpath2 = fullfile(x.fpath,username,fpath2);


%% Modify from here 
x.fname = 'test';
x.figure_on = 0;
x.file_type = 'TDT';

% Name of projects
x.projects = {
    'HP_timescale'
    };



% sampling rate. 

% blackrock or DBC
% x.fs = 30000;
% TDT
x.fs = 24414;




% Probe channel mapping information
x.Nb_ch = 32;
load('chanMap_TDT.mat');

x.chanMap = chanMap;
x.xcoords = round(xcoords/50);
x.ycoords = round(ycoords/50);


