for n = 1:length(GLM_dataset)
    for t = 1:length(GLM_dataset{n,3})
        if GLM_dataset{n,3}(t,2) == 1 && GLM_dataset{n,3}(t,3) ==1
            GLM_dataset{n,3}(t,6) = 1;
        elseif GLM_dataset{n,3}(t,2) == 1 && GLM_dataset{n,3}(t,3) ==0
            GLM_dataset{n,3}(t,6) = 0;
        elseif GLM_dataset{n,3}(t,2) == 2 && GLM_dataset{n,3}(t,3) ==0
            GLM_dataset{n,3}(t,6) = 1;
        elseif GLM_dataset{n,3}(t,2) == 2 && GLM_dataset{n,3}(t,3) ==1
            GLM_dataset{n,3}(t,6) = 0;
        end
    end
end