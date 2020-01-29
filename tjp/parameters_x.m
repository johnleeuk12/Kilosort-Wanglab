function x = parameters_x()


%Modify parameters based on OpenEphys and xbz output. 

addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
x.PC_name = '426_Analysis';
x.animal_name = 'M12E';
x.list  = {'2019-07-11_15-59-38'
'2019-07-11_16-01-27'
'2019-07-11_16-03-17'
'2019-07-11_16-05-08'
'2019-07-11_16-07-13'
'2019-07-11_16-09-23'
'2019-07-11_16-16-28'
'2019-07-11_16-25-00'
'2019-07-11_16-28-33'
'2019-07-11_16-36-06'
'2019-07-11_16-37-55'
'2019-07-11_16-40-23'
'2019-07-11_16-51-20'
'2019-07-11_16-59-55'
'2019-07-11_17-04-09'



};
x.fpath2 = 'D:\Data\Experiments\M12E\H3T1S1_concat';
x.session_id = 15;
x.xbz_file_name = 'M12E0289';


x.session_name = x.list{x.session_id};
x.file_type = '100';

x.fpath = directories(x.PC_name,x.animal_name,x.session_name);


x.fs = 30000;
x.Nb_ch = 64;
x.chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
    22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
    64 60 62 58 56 52 54 48 49 53 51 55 57 61 59 63 ...
    38 42 40 45 43 35 33 47 36 50 34 44 46 39 41 37];

