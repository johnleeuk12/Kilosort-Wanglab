% add to pool 
function create_GLM(Pool,name)


%{
UPDATE LOG
==========================================================================

20/01/2025 JHL
Adding transition phase to hit rate.
adding trialtype id to data


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


save(fullfile('D:\DATA',filesep,[name,'.mat']),"GLM_dataset")

% 
% sz = length(GLM_dataset);
% for p = 1:length(Pool)
%     GLM_dataset{p+sz,1} = Pool(p).spiketimes;
%     GLM_dataset{p+sz,2} = Pool(p).licktimes;
%     GLM_dataset{p+sz,4} = Pool(p).eventtimes;
%     GLM_dataset{p+sz,5} = Pool(p).best_ch;
%     idx = find(Pool(p).xb.trial_type(:,4) < 5);
%     M = zeros(length(idx),2);
%     for i = 1:length(idx)
%         if mod(Pool(p).xb.trial_type(i,4),2) ==1
%             M(i,1) = 1;
%             M(i,2) = 1;
%         elseif  mod(Pool(p).xb.trial_type(i,4),2) ==3
%             M(i,1) = 1;
%         elseif  mod(Pool(p).xb.trial_type(i,4),2) ==4
%             M(i,2) = 1;
%         end
%     end
%     
%         
%     GLM_dataset{p+sz,3} = [M(1:length(Pool(p).xb.hit_code),:),Pool(p).xb.hit_code];
% end