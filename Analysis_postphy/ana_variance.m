


gain ={};
stim = {};

clearvars -except Pool;


[gain{2},stim{2},rate,raster2] = ana_gain(Pool,2,1,0.5);
% [raster2, rate, Lick2, Lick_rate] = gather_raster_ephys(2, 6, Pool);

% for each unit in stim, calculate variance within R1, RT and R2
%%
r_ind = 2;
s_ind = 2; % 1: 5khz 2:10kHz

stim_list = find(stim{r_ind}{s_ind}.all ==1);
gain_list = find(gain{r_ind}{s_ind}.pv >0);
% stim_list = gain_list;
V = zeros(3,length(Pool));
M = zeros(3,length(Pool));

t1 = 2000;
t2 = 2500;
for n = stim_list
    if s_ind ==2
        V(1,n) = std(mean((rate.PSTH{n,2}(:,t1:t2)),2))^2;
        V(2,n) = std(mean((rate.PSTH{n,8}(:,t1:t2)),2))^2;
        V(3,n) = std(mean((rate.PSTH{n,4}(:,t1:t2)),2))^2;
        M(1,n) = mean2(rate.PSTH{n,2}(:,t1:t2));
        M(2,n) = mean2(rate.PSTH{n,8}(:,t1:t2));
        M(3,n) = mean2(rate.PSTH{n,4}(:,t1:t2));
    elseif s_ind ==1
        V(1,n) = std(mean((rate.PSTH{n,1}(:,t1:t2)),2))^2;
        V(2,n) = std(mean((rate.PSTH{n,7}(:,t1:t2)),2))^2;
        V(3,n) = std(mean((rate.PSTH{n,3}(:,t1:t2)),2))^2;
        M(1,n) = mean2(rate.PSTH{n,1}(:,t1:t2));
        M(2,n) = mean2(rate.PSTH{n,7}(:,t1:t2));
        M(3,n) = mean2(rate.PSTH{n,3}(:,t1:t2));
    end
end



% for ind = 1:3
%     figure(ind)
%     y1 = V(ind,stim_list);
%     x =M(ind,stim_list);
%     mdl =fitlm(x,y1,'Intercept',false);
%     scatter(x,y1);
%     hold on 
%     plot(mdl)
%     title(mdl.Coefficients.Estimate)
%     ylim([0,200])
%     xlim([0,50])
% 
% 
% end


% figure
% boxplot(V(2:3,stim_list).'-V(1,stim_list).')
% % ylim([0.5,2])

F = zeros(3,length(stim_list));
for ind = 1:3
    F(ind,:) = V(ind,stim_list)./M(ind,stim_list);
end
figure
boxplot(F.')
ylim([0.5,10])
% [T,p] = kstest(F(3,:)-F(2,:))
[p,T] = signrank(F(1,:),F(2,:))
[p,T] = signrank(F(3,:),F(2,:))
[p,T] = signrank(F(3,:),F(1,:))


% figure
% edges = 0:5:120;
% histogram(V(1,stim_list),edges)
% hold on
% histogram(V(2,stim_list),edges)
% histogram(V(3,stim_list),edges)


%%


r_ind = 2;
s_ind = 2; % 1: 5khz 2:10kHz

stim_list = find(stim{r_ind}{s_ind}.all ==1);
gain_list = find(gain{r_ind}{s_ind}.pv >0);
V = zeros(3,length(Pool));

S = {};

phases = [2,8,4];
for n = stim_list
    temp_S = {};
    for ind = 1:length(phases)
        S{n,ind} = [];
        for tr = 1:raster2{n}(phases(ind)).nreps
            temp_ISI = raster2{n}(phases(ind)).spikes(raster2{n}(phases(ind)).rep == tr);
            temp_ISI = temp_ISI(2:end)-temp_ISI(1:end-1);
            S{n,ind} = [S{n,ind},temp_ISI];
        end
    end
end

CV = zeros(3,length(Pool));    
for n = stim_list
    for ind = 1:length(phases)
        CV(ind,n) = std(S{n,ind})/mean(S{n,ind});
    end
end
    
CV = CV(:,stim_list);    
figure
boxplot(CV.')

[T,p] = ttest2(CV(1,:),CV(2,:))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

