function x = parameters_distance()


x.PC_name = '426_Analysis';
x.animal_name = 'M12E';
addpath(genpath(fullfile('D:\Data\M12E\Units')));
addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\Analysis')
x.PC_name = '426_Analysis';
x.animal_name = 'M12E';
x.fpath1 = directories(x.PC_name,x.animal_name,[]);



x.hole_number = 6;
x.track_number = 2;
x.depth = 2870;
x.hemi = 'L';
x.segment_list = {'M12E0463'
    'M12E0464'
    'M12E0465'
    'M12E0466'
    'M12E0467'
    'M12E0468'
    'M12E0469'
    'M12E0470'
    'M12E0471'
    'M12E0472'
    'M12E0473'
    'M12E0474'
    'M12E0475'
    'M12E0476'
    'M12E0477'
    };

