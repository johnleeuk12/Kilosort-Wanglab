function x = parameters_distanceS()


x.PC_name = '426_Analysis';
x.animal_name = 'M160E';
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
        
    case '426_Analysis_CI'
        
end



x.hole_number = 2;
x.track_number = 1;
x.depth = 905;
x.hemi = 'L';
x.segment_list = {    
'M160E0024'
'M160E0025'
'M160E0026'
'M160E0027'
'M160E0028'
'M160E0029'
'M160E0030'

};







