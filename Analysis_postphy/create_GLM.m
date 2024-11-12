% add to pool 
function GLM_dataset = create_GLM(GLM_dataset,Pool)




sz = length(GLM_dataset);
for p = 1:length(Pool)
    GLM_dataset{p+sz,1} = Pool(p).spiketimes;
    GLM_dataset{p+sz,2} = Pool(p).licktimes;
    GLM_dataset{p+sz,4} = Pool(p).eventtimes;
    GLM_dataset{p+sz,5} = Pool(p).best_ch;
    idx = find(Pool(p).xb.trial_type(:,4) < 5);
    M = zeros(length(idx),2);
    for i = 1:length(idx)
        if mod(Pool(p).xb.trial_type(i,4),2) ==1
            M(i,1) = 1;
            M(i,2) = 1;
        elseif  mod(Pool(p).xb.trial_type(i,4),2) ==3
            M(i,1) = 1;
        elseif  mod(Pool(p).xb.trial_type(i,4),2) ==4
            M(i,2) = 1;
        end
    end
    
        
    GLM_dataset{p+sz,3} = [M(1:length(Pool(p).xb.hit_code),:),Pool(p).xb.hit_code];
end