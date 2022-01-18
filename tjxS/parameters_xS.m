function x = parameters_xS(PC, animal, path2, SID)

%Modify parameters based on OpenEphys and xbz output. 
addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')

x.PC_name = PC;

x.animal_name = animal;

x.list  = {
'2021-12-19_14-52-43'
'2021-12-19_14-59-23'
'2021-12-19_15-06-01'
'2021-12-19_15-12-43'
'2021-12-19_15-18-59'
'2021-12-19_15-27-48'
'2021-12-19_16-10-42'
'2021-12-19_16-25-34'
'2021-12-19_16-44-03'
'2021-12-19_17-00-38'
'2021-12-19_17-06-38'


};

x.fpath2 = path2;
x.session_id = SID;
x.figure_on = 0;
xbz_list = {
'M56E0052'
'M56E0053'
'M56E0054'
'M56E0055'
'M56E0056'
'M56E0057'
'M56E0059'
'M56E0060'
'M56E0061'
'M56E0062'
'M56E0063'

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

