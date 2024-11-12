function x = parameters_xS(PC, username,fpath2)

%% Do not change
x.PC_name = PC;
x.animal_name = username;
[x.fpath, x.savepath] = directories(x.PC_name);
x.fpath2 = fullfile(x.fpath,username,fpath2);


%% Modify from here 
x.fname = 'continuous';
x.figure_on = 0;
x.file_type = 'DBC';

% Name of projects
x.projects = {
    'HP_timescale', ...
    'RL'
    };
%     



% sampling rate. 

% blackrock or DBC
x.fs = 30000;
% TDT
% x.fs = 24414;




% Probe channel mapping information
x.Nb_ch = 64;
x.Nb_ch_real = 64;
% load('
load('chanMap_DBC_NN.mat');
% 
% x.chanMap = chanMap(1:x.Nb_ch_real);
% x.xcoords = round(xcoords(1:x.Nb_ch_real)/50);
% x.ycoords = round(ycoords(1:x.Nb_ch_real)/50);


x.chanMap = chanMap;
x.xcoords = round(xcoords/50);
x.ycoords = round(ycoords/50);

