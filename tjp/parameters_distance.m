function x = parameters_distance()


x.PC_name = '426_John';
% x.PC_name = '426_Analysis';
x.animal_name = 'M12E';
addpath(genpath(fullfile('D:\Data\M12E\Units')));
% addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
% x.animal_name = 'M12E';
% x.fpath1 = directories(x.PC_name,x.animal_name,[]);



x.hole_number = 13;
x.track_number = 6;
x.depth = 800;
x.hemi = 'L';
x.segment_list = {
    
'M12E1417'
'M12E1418'
'M12E1419'
'M12E1420'
'M12E1421'
'M12E1422'


    };

