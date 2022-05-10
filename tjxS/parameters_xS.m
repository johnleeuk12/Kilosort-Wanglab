function x = parameters_xS(PC, animal, path2, SID)

%Modify parameters based on OpenEphys and xbz output. 
addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')

x.PC_name = PC;

x.animal_name = animal;

x.list  = {
'2022-02-15_15-09-11'
'2022-02-15_15-18-00'
'2022-02-15_15-26-19'
'2022-02-15_15-34-45'
'2022-02-15_15-48-41'
'2022-02-15_15-57-20'
'2022-02-15_16-12-08'


};


x.fpath2 = path2;
x.session_id = SID;
x.figure_on = 0;
xbz_list = {
'M56E0520'
'M56E0521'
'M56E0522'
'M56E0523'
'M56E0524'
'M56E0525'
'M56E0526'

};


x.xbz_file_name = xbz_list{x.session_id};

x.session_name = x.list{x.session_id};

x.file_type = '100';

[x.fpath, x.savepath] = directories(x.PC_name,x.animal_name,x.session_name);

x.fs = 30000;
x.Nb_ch = 64;
x.chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
    22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
    64 59 61 57 55 51 53 49 48 54 52 56 58 62 60 64 ...
    37 41 39 46 44 34 50 36 47 33 35 43 45 40 42 38];

