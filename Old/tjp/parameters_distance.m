function x = parameters_distance()


x.PC_name = '426_John';
% x.PC_name = '426_Analysis';
x.animal_name = 'M60F';
addpath(genpath(fullfile('D:\Data\Units\M60F')));
% addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
% x.animal_name = 'M12E';
% x.fpath1 = directories(x.PC_name,x.animal_name,[]);



x.hole_number = 14;
x.track_number =6;
x.depth = 800;
x.hemi = 'L';
x.segment_list = {
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



    