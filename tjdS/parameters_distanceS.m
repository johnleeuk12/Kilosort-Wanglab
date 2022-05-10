function x = parameters_distanceS()


x.PC_name = '426_Analysis';
x.animal_name = 'M56E';
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



x.hole_number = 5;
x.track_number = 8;
x.depth = 905;
x.hemi = 'L';
x.segment_list = { 
'M56E0520'
'M56E0521'
'M56E0522'
'M56E0523'
'M56E0524'
'M56E0525'
'M56E0526'

};







