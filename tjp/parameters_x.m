function x = parameters_x()


%Modify parameters based on OpenEphys and xbz output. 

addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
x.PC_name = '426_Analysis';
x.animal_name = 'M12E';
x.list  = {'2019-07-26_16-31-51'
'2019-07-26_16-34-50'
'2019-07-26_16-42-16'
'2019-07-26_16-47-30'
'2019-07-26_16-56-03'
'2019-07-26_17-05-34'
'2019-07-26_17-08-48'
'2019-07-26_17-15-29'
'2019-07-26_17-21-53'
'2019-07-26_17-26-29'
'2019-07-26_17-31-04'
'2019-07-26_17-32-54'
'2019-07-26_17-34-41'
'2019-07-26_17-36-30'
'2019-07-26_17-38-31'
};

x.session_id = 1;
x.session_name = x.list{x.session_id};
x.file_type = '100';
x.xbz_file_name = 'M12E0463';
x.fpath = directories(x.PC_name,x.animal_name,x.session_name);

x.fpath2 = 'C:\Data\OpenEphys\M12E\H6T3S1_concat';
x.fs = 30000;
x.Nb_ch = 64;
x.chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
    22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
    64 60 62 58 56 52 54 48 49 53 51 55 57 61 59 63 ...
    38 42 40 45 43 35 33 47 36 50 34 44 46 39 41 37];

