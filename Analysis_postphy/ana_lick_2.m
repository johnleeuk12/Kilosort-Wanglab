
load('D:\DATA\Pool_PPC4.mat');
% [raster2, rate, Lick2, Lick_rate] = gather_raster_ephys(2, 6, Pool,list_lick.');

% list_lick = list_lick+1;

%%

r_data ={};
delta = zeros(length(Pool),2);

soft = 0.5;
for n = 1:length(list_lick)
    nn = list_lick(n);
    hit_code = Pool(nn).xb.hit_code;
    if length(hit_code)>460
        hit_code = hit_code(1:460,:);
    end
    r_data{1} = find(hit_code(:,1) == 1); % hit 
    r_data{1} = cat(2,r_data{1},zeros(length(r_data{1}),5));
    
    r_data{2} = find(hit_code(:,4) == 1); % FA
    r_data{2} = cat(2,r_data{2},zeros(length(r_data{2}),5));
    for ind = 1:2
        for t = 1:length(r_data{ind}(:,1))
            tr = r_data{ind}(t,1);
            st_onset = Pool(nn).xb.trial_type(tr,1)*1e-3;
            lick = Pool(nn).licktimes(find(Pool(nn).licktimes>st_onset));
            r_onset = lick((lick>st_onset + 1.5));
            r_data{ind}(t,2)= min(r_onset);
            r_data{ind}(t,3) = sum((Pool(nn).spiketimes < r_data{ind}(t,2) & Pool(nn).spiketimes > r_data{ind}(t,2) -0.25));
            r_data{ind}(t,4) = sum((Pool(nn).spiketimes < r_data{ind}(t,2)+0.35 & Pool(nn).spiketimes > r_data{ind}(t,2) +0.10));
        end
        
        delta(nn,ind) = (mean(r_data{ind}(:,4))- mean(r_data{ind}(:,3)))/...
            (max([mean(r_data{ind}(:,4)),mean(r_data{ind}(:,3))])+soft);
    end
    
    

end
    