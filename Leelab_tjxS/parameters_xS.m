function x = parameters_xS(PC, animal,fpath2)

%Modify parameters based on OpenEphys and xbz output. 
addpath('D:\GitHub\Kilosort-Wanglab\Analysis')


addpath(genpath('D:\GitHub\npy-matlab'))
addpath(genpath('D:\GitHub\TDTMatlabSDK'))
x.PC_name = PC;
x.animal_name = animal;


x.fname = 'TDTdata.dat';
x.figure_on = 0;
x.file_type = 'TDT';

[x.fpath, x.savepath] = directories(x.PC_name,x.animal_name);

x.fpath2 = fpath2;

% blackrock or DBC
% x.fs = 30000;



% TDT
x.fs = 24414;
x.Nb_ch = 32;
load('chanMap_TDT.mat');

x.chanMap = chanMap;
x.xcoords = round(xcoords/50);
x.ycoords = round(ycoords/50);
% x.chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
%     22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
%     64 60 62 58 56 52 54 48 49 53 51 55 57 61 59 63 ...
%     38 42 40 45 43 35 33 47 36 50 34 44 46 39 41 37];

