function x = parameters_xS(PC, animal, path2, SID)

%Modify parameters based on OpenEphys and xbz output. 
addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')

x.PC_name = PC;

x.animal_name = animal;

x.list  = {
'2021-01-02_15-31-30'
'2021-01-02_15-47-34'
'2021-01-02_16-02-19'
'2021-01-02_16-22-30'
'2021-01-02_16-40-31'
'2021-01-02_16-53-55'
'2021-01-02_16-58-46'
'2021-01-02_17-03-04'
'2021-01-02_17-07-22'
'2021-01-02_17-11-39'
'2021-01-02_17-15-21'

};

x.fpath2 = path2;
x.session_id = SID;
x.figure_on = 0;
xbz_list = {
'M60F0273'
'M60F0274'
'M60F0275'
'M60F0276'
'M60F0277'
'M60F0278'
'M60F0279'
'M60F0280'
'M60F0281'
'M60F0282'
'M60F0283'

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

