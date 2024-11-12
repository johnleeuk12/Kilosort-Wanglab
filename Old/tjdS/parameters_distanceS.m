function x = parameters_distanceS()


x.PC_name = 'Leelab_John';
x.animal_name = 'M60F';
% x.list_data_name = [x.animal_name '_unit_list']; %Do not modify this
x.list_data_name = 'unit_list'; %Do not modify this


switch x.PC_name
    case '426_Analysis'
        addpath('C:\Users\Seth\Documents\GitHub\Kilosort-Wanglab\tjdS');
        x.save_dir = fullfile('D:\Data\Units', filesep, x.animal_name);
        try
            addpath(x.save_dir);
        catch
            error('path not found')
        end
        
    case 'Leelab_John'
        addpath('D:\GitHub\Kilosort-Wanglab\tjdS');
        x.save_dir = fullfile('D:\Data\Units', filesep, x.animal_name);
        try
            addpath(x.save_dir);
        catch
            error('path not found')
        end
        
end



x.hole_number = 3;
x.track_number = 6;
x.depth = 2180;
x.hemi = 'L';
x.segment_list = {    
'M60F0161'
'M60F0162'
'M60F0163'
'M60F0164'

};







