function x = parameters_x()


%Modify parameters based on OpenEphys and xbz output. 

addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
x.PC_name = '426_Analysis';
x.animal_name = 'M12E';
x.list  = {'2019-07-27_14-24-54'
'2019-07-27_14-26-39'
'2019-07-27_14-28-36'
'2019-07-27_14-30-40'
'2019-07-27_14-33-57'
'2019-07-27_14-38-11'
'2019-07-27_14-40-28'
'2019-07-27_14-50-19'
'2019-07-27_14-59-45'
'2019-07-27_15-07-30'
'2019-07-27_15-10-49'
'2019-07-27_15-18-16'
'2019-07-27_15-21-46'
'2019-07-27_15-26-37'
'2019-07-27_15-30-37'




};

x.session_id = 15;
x.session_name = x.list{x.session_id};
x.file_type = '100';
x.xbz_file_name = 'M12E0492';
x.fpath = directories(x.PC_name,x.animal_name,x.session_name);

x.fpath2 = 'D:\Data\Experiments\M12E\H6T4S1_concat';
x.fs = 30000;
x.Nb_ch = 64;
x.chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
    22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
    64 60 62 58 56 52 54 48 49 53 51 55 57 61 59 63 ...
    38 42 40 45 43 35 33 47 36 50 34 44 46 39 41 37];
