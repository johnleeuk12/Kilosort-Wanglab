function x = parameters_xS(PC, animal, path2, SID)

%Modify parameters based on OpenEphys and xbz output. 
addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')

x.PC_name = PC;

x.animal_name = animal;

x.list  = {
'2021-05-10_15-14-10'
'2021-05-10_15-19-36'
'2021-05-10_15-24-54'
'2021-05-10_15-30-10'
'2021-05-10_15-35-29'
'2021-05-10_15-51-49'
'2021-05-10_16-07-42'
'2021-05-10_16-25-06'
'2021-05-10_16-39-45'

};



x.fpath2 = path2;
x.session_id = SID;
x.figure_on = 0;
xbz_list = {
'M60F1141'
'M60F1142'
'M60F1143'
'M60F1144'
'M60F1145'
'M60F1146'
'M60F1147'
'M60F1148'
'M60F1149'

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

