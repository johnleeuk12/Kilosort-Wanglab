addpath('D:\Data\M12E\Units_training')

x = dir;

for i = 1:length(x)
    if x(i).bytes ~= 0 && ~strcmp(x(i).name,'M12E_unit_list.mat')
        %         disp(i)
        
        unit_name = x(i).name;
        %         if ~strcmp(unit_name(end-4),'v')
        %             if strcmp(unit_name(5),'u')
        %                 unit_name = ['M12Eu00' unit_name(end-6:end)];
        %             else
        %                 unit_name = ['M12Eu0' unit_name(end-7:end)];
        %             end
        %
        %         end
        unit_n = num2str(str2num(unit_name(6:end-4))+9641);
        if length(unit_n) == 4
            unit_name = ['M12Eu0',unit_n,'.mat'];
        else
            unit_name = ['M12Eu',unit_n,'.mat'];
        end
        dos(['rename "' x(i).name '" "' unit_name '"']); % (1)
    end
end