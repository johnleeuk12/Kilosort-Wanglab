function x = parameters_x()

%Modify parameters based on OpenEphys and xbz output. 

addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
x.PC_name = '426_Analysis';
x.animal_name = 'M12E';
x.list  = {
'2020-01-04_16-08-52'
'2020-01-04_16-12-52'
'2020-01-04_16-15-23'
'2020-01-04_16-19-42'
'2020-01-04_16-22-15'
'2020-01-04_16-25-41'
'2020-01-04_16-36-01'
'2020-01-04_16-37-53'
'2020-01-04_16-40-06'

};
x.fpath2 = 'D:\Data\Experiments\M12E\H9T2S1_concat';
x.session_id = 1;
x.figure_on = 1;
xbz_list = {
'M12E0842'
'M12E0843'
'M12E0844'
'M12E0845'
'M12E0846'
'M12E0847'
'M12E0848'
'M12E0849'
'M12E0850'

};



x.xbz_file_name = xbz_list{x.session_id};

x.session_name = x.list{x.session_id};

x.file_type = '100';

x.fpath = directories(x.PC_name,x.animal_name,x.session_name);


x.fs = 30000;
x.Nb_ch = 64;
x.chanMap = [27 32 21 3 25 30 19 5 23 17 24 7 20 29 26 9 ...
    22 31 28 11 16 13 1 15 18 12 14 10 8 4 6 2 ...
    64 60 62 58 56 52 54 48 49 53 51 55 57 61 59 63 ...
    38 42 40 45 43 35 33 47 36 50 34 44 46 39 41 37];

