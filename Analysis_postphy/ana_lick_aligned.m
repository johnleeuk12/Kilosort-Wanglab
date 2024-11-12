function [C2, spont, L] = ana_lick_aligned(Pool,rule,ops)



%% 06242024 Align to nth lick

% bin = 5;
% edges = -2e3:bin:4e3;

bin = ops.bin;
edges = ops.edges;
pre = 2;
post = 6;


C2 = {};
C2.hit = zeros(length(Pool),length(edges)-1);
C2.hit =C2.hit.';
C2.FA = zeros(length(Pool),length(edges)-1);
C2.FA =C2.FA.';
C2.all = zeros(length(Pool),length(edges)-1);
C2.all = C2.all.';


L = {};
L.hit = zeros(length(edges)-1,length(Pool));
L.FA = zeros(length(edges)-1,length(Pool));

spont = {};
spont.mean = zeros(length(Pool),1);
spont.hit = zeros(length(Pool),1);
spont.FA = zeros(length(Pool),1);
spont.std = zeros(length(Pool),3);


tic
for n = 1:length(Pool)
    if mod(n,50) ==1
        fprintf(['%4d /' num2str(length(Pool)) ' time : %6.2f sec \n'],n,toc')
    end
    
    try
        event_times= Pool(n).eventtimes(Pool(n).xb.trial_type(1:length(Pool(n).eventtimes),4) < 5);
    catch
        event_times= Pool(n).eventtimes(Pool(n).xb.trial_type(:,4)< 5);
    end

    ILI = Pool(n).licktimes(2:end)-Pool(n).licktimes(1:end-1);
    lick_onset= find(ILI>2); % lick bout = 2 seconds
    lick_onset = lick_onset+1;
    lick_offset = lick_onset-1;
    
    % event_times = event_times + 1.5; % add 1.5 seconds to eventtimes to align to reward period.
    Rt = zeros(length(Pool(n).xb.hit_code),10);
    spont_raw = zeros(length(event_times),1);
    Spks = {};
    Lcks = {};
    for tr  = 1:length(event_times)
        lick_ind =  find(Pool(n).licktimes > event_times(tr) & ...
                        Pool(n).licktimes < (event_times(tr) + post) ); 
        
        rt_ind = find(Pool(n).licktimes(lick_onset) > event_times(tr) & ...
                        Pool(n).licktimes(lick_onset) < (event_times(tr) + post) );
        sp_ind = find(Pool(n).spiketimes > event_times(tr) & ...
                        Pool(n).spiketimes < (event_times(tr)+post) );
                    
        sp_ind2 = find(Pool(n).spiketimes > event_times(tr)-2 & ...
                        Pool(n).spiketimes < (event_times(tr)-0.5) ); % spont rate, 1 sec before event
        if length(rt_ind) >10
            Rt(tr,:) = rt_ind(1:10);
        else
            Rt(tr,1:length(rt_ind)) = rt_ind;
        end
        Spks{tr,1} = Pool(n).spiketimes(sp_ind);% -event_times(tr);
        spont_raw(tr,1) = length(Pool(n).spiketimes(sp_ind2));
        Lcks{tr,1} = Pool(n).licktimes(lick_ind);
    end
    
    FR = {};

    FR.Hit = [];
    FR.FA = [];
    FR.all = [];
    LR = {};
    LR.Hit = [];
    LR.FA = [];
    LR.all = [];
    
    switch rule
        case 1
            hit_ind = find(Pool(n).xb.hit_code(1:200,1)==1);
            FA_ind = find(Pool(n).xb.hit_code(1:200,3)==1);
        case 2
            hit_ind = find(Pool(n).xb.hit_code(201:end,1)==1);
            FA_ind = find(Pool(n).xb.hit_code(201:end,3)==1);
            hit_ind = hit_ind +200;
            FA_ind = FA_ind+200;
        case 0
            hit_ind = find(Pool(n).xb.hit_code(:,1)==1);
            FA_ind = find(Pool(n).xb.hit_code(:,3)==1);
            
    end
    



  
    for r = 1:length(hit_ind)
        for ln = 1:sum(Rt(hit_ind(r),:)>0)
            LR.Hit = [LR.Hit; Lcks{hit_ind(r),1}-Pool(n).licktimes(lick_onset(Rt(hit_ind(r),1)))];
            try
                FR.Hit = [FR.Hit;Spks{hit_ind(r)}-Pool(n).licktimes(lick_onset(Rt(hit_ind(r),ln)))];
            catch
            end
        end
    end
    for r = 1:length(FA_ind)
        for ln = 1:sum(Rt(FA_ind(r),:)>0)
            LR.FA = [LR.FA; Lcks{FA_ind(r),1}-Pool(n).licktimes(lick_onset(Rt(FA_ind(r),1)))];           
            try
                FR.FA = [FR.FA;Spks{FA_ind(r)}-Pool(n).licktimes(lick_onset(Rt(FA_ind(r),ln)))];
            catch
            end
        end
        
    end
    for tr  = 1:length(event_times)
        for ln = 1:sum(Rt(tr,:)>0)
            LR.all = [LR.all; Lcks{tr,1}-Pool(n).licktimes(lick_onset(Rt(tr,1)))];
            try
                FR.all = [FR.Hit;Spks{tr}-Pool(n).licktimes(lick_onset(Rt(tr,ln)))];
            catch
            end
        end
    end
    
    nb_hit2 = sum(sum(Rt(hit_ind,:)>0)); 
    nb_FA2 = sum(sum(Rt(FA_ind,:)>0));
    

    countM = histcounts(FR.Hit*1e3,edges);
    C2.hit(:,n) = countM*1e3/(bin*nb_hit2);
    countM = histcounts(FR.FA*1e3,edges);
    C2.FA(:,n) = countM*1e3/(bin*nb_FA2);
    countM = histcounts(FR.Hit*1e3,edges);
    C2.all(:,n) = countM*1e3/(bin*nb_hit2);
    spont.mean(n,1) = mean(spont_raw); % in FR, since we're counting nb of spks per sec
%     spont.mean(n,2) = mean(spont_raw(hit_ind));
    spont.std(n,1) = std(spont_raw,[],1);
%     spont.std(n,2) = std(spont_raw(hit_ind),[],1);

    L.hit(:,n) = histcounts(LR.Hit*1e3,edges)*1e3/(bin*nb_hit2);
    L.FA(:,n) = histcounts(LR.FA*1e3,edges)*1e3/(bin*nb_FA2);
  
end