% add to pool 
function create_GLM(Pool,name,D)


%{
UPDATE LOG
==========================================================================

20/01/2025 JHL
Adding transition phase to hit rate.
adding trialtype id to data

14/07/2025 JHL

before running the code, load D as the following: 
D = load('D:\DATA\List data\ICPool_Event_TTR.mat');


%}

% sz = length(GLM_dataset);
for p = 1:length(Pool)
    GLM_dataset{p,1} = Pool(p).spiketimes;
    GLM_dataset{p,2} = Pool(p).licktimes;
    GLM_dataset{p,4} = Pool(p).xb.hit_code;
    GLM_dataset{p,3} = Pool(p).eventtimes;
    GLM_dataset{p,5} = Pool(p).xb.trial_type;
    GLM_dataset{p,6} = Pool(p).best_ch;

end

%% Adding TTR and Lickbout onset times

dataset= D.dataset;

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

%%

save(fullfile('D:\DATA',filesep,[name,'.mat']),"GLM_dataset")

