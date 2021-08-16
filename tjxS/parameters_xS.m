function x = parameters_xS(PC, animal, path2, SID)

%Modify parameters based on OpenEphys and xbz output. 
addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')

x.PC_name = PC;

x.animal_name = animal;

x.list  = {
'2021-08-07_15-34-10'
'2021-08-07_15-40-49'
'2021-08-07_15-47-34'
'2021-08-07_15-56-43'
'2021-08-07_16-17-10'
'2021-08-07_16-33-36'
'2021-08-07_16-48-48'
'2021-08-07_17-03-57'
'2021-08-07_17-10-57'
'2021-08-07_17-17-42'



};

x.fpath2 = path2;
x.session_id = SID;
x.figure_on = 0;
xbz_list = {
'M160E0036'
'M160E0037'
'M160E0038'
'M160E0039'
'M160E0040'
'M160E0041'
'M160E0042'
'M160E0043'
'M160E0044'
'M160E0045'

};

x.xbz_file_name = xbz_list{x.session_id};

x.session_name = x.list{x.session_id};

x.file_type = '100';

[x.fpath, x.savepath] = directories(x.PC_name,x.animal_name,x.session_name);

x.fs = 30000;
x.Nb_ch = 64;
x.chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
    22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
    64 60 62 58 56 52 54 48 49 53 51 55 57 61 59 63 ...
    38 42 40 45 43 35 33 47 36 50 34 44 46 39 41 37];

