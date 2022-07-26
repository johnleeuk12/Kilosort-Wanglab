

addpath('D:\GitHub\npy-matlab')
addpath('D:\GitHub\Kilosort-Wanglab\Analysis')
x.fpath2 = 'D:\Data\example_TDT recording data\';
x.fname = 'TDTdata.dat';

load(fullfile(x.fpath2,'chanmap.mat'));
x.chanMap = chanMap;
x.fs = 24414;
x.Nb_ch = 32;
obj.params = x;