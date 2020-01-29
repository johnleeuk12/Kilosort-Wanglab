function x = parameters_distance()


x.PC_name = '426_John';
% x.PC_name = '426_Analysis';
x.animal_name = 'M12E';
addpath(genpath(fullfile('D:\Data\M12E\Units')));
% addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
% x.animal_name = 'M12E';
% x.fpath1 = directories(x.PC_name,x.animal_name,[]);



x.hole_number = 2;
x.track_number = 7;
x.depth = 2900;
x.hemi = 'L';
x.segment_list = {'M12E0260'
'M12E0261'
'M12E0262'
'M12E0263'
'M12E0264'
'M12E0265'
'M12E0266'
'M12E0267'
'M12E0268'


    };

