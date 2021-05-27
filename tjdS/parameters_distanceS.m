function x = parameters_distanceS()


x.PC_name = '426_Analysis';
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
        
    case '426_Analysis_CI'
        
end



x.hole_number = 15;
x.track_number = 2;
x.depth = 1300;
x.hemi = 'L';
x.segment_list = {    
'M60F1194'
'M60F1195'
'M60F1196'
'M60F1197'
'M60F1198'
'M60F1199'
'M60F1200'
'M60F1201'



};







