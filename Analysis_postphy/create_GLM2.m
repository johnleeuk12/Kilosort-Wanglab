%% Create GLM 2

%{
Function to integrate lickbout onset times to trial type.
using Pool_data_TTR

%}

for n = 1:length(GLM_dataset)
    lick_bouts = dataset{n,2}((dataset{n,2}(:,3)>0),1)*1e-3;
    for tr = 1:length(GLM_dataset{n,3})
        if tr == length(GLM_dataset{n,3})
            diff = lick_bouts(find(lick_bouts>GLM_dataset{n,3}(tr,1)));
        else
            diff = lick_bouts(find(lick_bouts>GLM_dataset{n,3}(tr,1) & lick_bouts<GLM_dataset{n,3}(tr+1,1)));
        end
        if length(diff)>1
            diff = diff(1);
        end
        if isempty(diff)
            GLM_dataset{n,3}(tr,2) = 0;
        else
            GLM_dataset{n,3}(tr,2) = diff-GLM_dataset{n,3}(tr,1);
        end
        if  GLM_dataset{n,3}(tr,2) > 8
            GLM_dataset{n,3}(tr,2) = 0;
        end
        
    end
    
    GLM_dataset{n,7} = dataset{n,5}+200;
end

%% TTR

for n = 1:length(GLM_dataset)
    ttr = dataset{n,3}(:,1)- dataset{n,4}(1,1);
    ttr = abs(ttr);
    [~,ind] = min(ttr);
    GLM_dataset{n,7} = ind+1;
end

%% add TTR to pool


D = load('D:\DATA\List data\ICPool_Event_TTR.mat');
dataset= D.dataset;

for n = 1:length(Pool)
    ttr = dataset{n,3}(:,1)- dataset{n,4}(1,1);
    ttr = abs(ttr);
    [~,ind] = min(ttr);
    Pool(n).ttr = ind+1;
end


%% add waveforms to Pool























