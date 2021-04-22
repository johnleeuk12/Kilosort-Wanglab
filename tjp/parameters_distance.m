function x = parameters_distance()


x.PC_name = '426_John';
% x.PC_name = '426_Analysis';
x.animal_name = 'M60F';
addpath(genpath(fullfile('D:\Data\Units\M60F')));
% addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
% x.animal_name = 'M12E';
% x.fpath1 = directories(x.PC_name,x.animal_name,[]);



x.hole_number = 5;
x.track_number = 4;
x.depth = 800;
x.hemi = 'L';
x.segment_list = {
'M60F0388'
'M60F0389'
'M60F0390'
'M60F0391'


    };

