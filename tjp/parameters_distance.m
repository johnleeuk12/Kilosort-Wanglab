function x = parameters_distance()


x.PC_name = '426_John';
% x.PC_name = '426_Analysis';
x.animal_name = 'M60F';
addpath(genpath(fullfile('D:\Data\Units\M60F')));
% addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
% x.animal_name = 'M12E';
% x.fpath1 = directories(x.PC_name,x.animal_name,[]);



x.hole_number = 6;
x.track_number = 7;
x.depth = 800;
x.hemi = 'L';
x.segment_list = {
'M60F0467'
'M60F0468'
'M60F0469'
'M60F0470'
'M60F0471'
'M60F0472'
'M60F0473'
'M60F0474'
'M60F0475'

    };

