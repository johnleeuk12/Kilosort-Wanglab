% Analysis of stim encoding 

% plot(smooth(mean(rate.PSTH{96,4},1),10))

n = 118;
W = zeros(4,length(rate.PSTH));

for n = 1:length(Pool)
    R = [rate.PSTH{n,1};rate.PSTH{n,2};rate.PSTH{n,3};rate.PSTH{n,4}];

    R2 = zeros(4,length(R));
%     figure
    for f = 1:4
        R2(f,:) = mean(rate.PSTH{n,f},1)-mean(R,1);
        R2(f,:) = smooth(R2(f,:),10);

%         plot(smooth(mean(R2(f,:),1),100))
%         hold on
    end
    R2 = R2/(max(max(abs(R2)))+0.5);
% 
    figure
    for f = 1:4
        plot(smooth(mean(R2(f,:),1),1))
        hold on
    end
%     % R2 = rate.PSTH{n,4}-mean(R,1);
%     % R2 = R2/max((abs(smooth(mean(R2,1),10))))+0.5;

    for f = 1:4
        W(f,n) = mean(R2(f,2000:2200));
    end
end
% plot(smooth(mean(R2,1),100))

figure
scatter(W(1,:),W(3,:))
hold on
xline([0,0.1,-0.1])
yline([0,0.1,-0.1])
xlim([-1,1]);
ylim([-1,1]);


% coeff = pca(W);

%% Calculate cross correlation 

% load('Pool_IC.mat');
% n = 155;
fig_on = 0;
L = zeros(1,length(Pool));
edges = [-2e2:5:3e2];
C = zeros(length(edges)-1,length(Pool));
tic

n_list = (1:length(Pool));
n_list = n_list(list3);
for n = n_list %1:length(Pool)
    if mod(n,10) ==1
        fprintf(['%4d /' num2str(length(Pool)) ' time : %6.2f sec \n'],n,toc')
    end
    
    M = zeros(length(Pool(n).spiketimes),length(Pool(n).licktimes));
    
    for s  = 1:length(Pool(n).spiketimes)
        M(s,:) = floor(Pool(n).spiketimes(s)*1e3)-floor(Pool(n).licktimes*1e3);
    end
    
    M = reshape(M,[],1);
    
    %     histogram(M,edges)
    countM = histcounts(M,edges);
    C(:,n) = countM;
    %     figure(n)
    %     plot(edges(1:end-1),countM);
    
    
    [pks,loc] = max(countM);
    if edges(loc) < 250 && loc >1
        pks2 = countM(loc-1:loc+1);
        p = prctile(countM,50); % 25th percentile
        sd= std(countM);
        
        if fig_on ==1
            figure(n)
            histogram(M,edges)
            hold on
            yline(2*sd + p,'--');
        end
        if sum((pks2 > 2*sd + p))>1
            
            L(1,n) = edges(loc);
            if fig_on == 1
                scatter(edges(loc+1),pks,50,'vr','filled')
            end
            
        end
    end
    %     plot(edges(1:end-1),countM)
    
    if fig_on ==1
        drawnow()
        pause(0.1)
    end

    
    
    
    
    
end


% edges = -2e2:20:2e2;
% histogram(L(L~=0),edges)
% ylim([0,10]);

%% with prelick list
C_norm = zeros(size(C,1),size(C,2));
for n = n_list
    C_norm(:,n) = (C(:,n)-prctile(C(:,n),25))/(max(C(:,n)-prctile(C(:,n),25))+10);
end

figure
plot(edges(1:end-1),mean(C_norm(:,n_list),2))

%% 

L = L_all.AC;
C = C_all.AC;
T = find(L ~= 0);
C_norm = zeros(size(C,1),size(C,2));
for n  = 1:length(C)
    C_norm(:,n) = (C(:,n)-prctile(C(:,n),25))/(max(C(:,n)-prctile(C(:,n),25))+10);
end

%%
figure
% plot(edges(1:end-1),mean(C_norm(:,T),2))
% C_norm_all{2} = C_norm(:,T)
c_list = {'r','g','b'};
for i = 1:3
    shadedErrorBar(edges(1:end-1),mean(C_norm_all{i},2),...
        std(C_norm_all{i},0,2)/sqrt(size(C_norm_all{i},2)),'lineProps',c_list{i})
    hold on
end

xlim([-200,250])


%% Separate Reward lick vs Spont lick


T = find(L ~= 0);
%%
edges = [-2e2:5:3e2];
C = {};
for ss = 1:8
    C{ss} = zeros(length(edges)-1,length(T));
end
for n = T
    M = {};
    for ss = 1:4
        M{ss} = [];
        for r = 1:80 % nreps, take first 80 reps
            sp_list = raster2{1,n}(ss).spikes(find(raster2{1,n}(ss).rep ==r));
            l_list = Lick2{1,n}(ss).licks(find(Lick2{1,n}(ss).rep ==r));
            Mp = zeros(length(sp_list),length(l_list));
            for s = 1:length(sp_list)
                Mp(s,:) = l_list - sp_list(s);
            end
            Mp = reshape(Mp,[],1);
            M{ss} = [M{ss};Mp];
            
        end
        M{ss} = M{ss}*1e3;
        countM = histcounts(M{ss},edges);
        C{ss}(:,n) = countM;
    
    end
    C{1}(:,n) = histcounts([M{1};M{4}],edges);
    C{2}(:,n) = histcounts([M{2};M{3}],edges);
%     C{1}(:,n) = histcounts([M{7}],edges);
%     C{2}(:,n) = histcounts([M{8}],edges);
% 
%     figure(n)
%     for ss = 1:4
%         subplot(2,2,ss)
%         histogram(M{ss},edges)
%     end
% 
%     pause(0.1)
%     drawnow()
end
    
        
figure
C_norm = {};
for ss = 1:2
    C_norm{ss} = zeros(size(C{ss},1),size(C{ss},2));
    for n  = 1:length(C{ss})
        C_norm{ss}(:,n) = (C{ss}(:,n)-prctile(C{ss}(:,n),25))/(max(C{ss}(:,n)-prctile(C{ss}(:,n),25))+5);
    end
    
    shadedErrorBar(edges(1:end-1),mean(C_norm{ss}(:,T),2),...
        std(C_norm{ss}(:,T),[],2)/sqrt(size(C_norm{ss}(:,T),2)),'lineProps',c_list{ss})
    hold on
end
        
for ss = 1:2
    figure
    s_list = zeros(1,length(T));
    if ss ==1
        for t = 1:length(T)
            [~,s_list(t)] = max(C_norm{ss}(:,T(t)));
        end
    end
    [~,I] = sort(s_list);
%     argmax(C_norm{ss}(:,T),1);
    imagesc(edges(1:end-1),1:length(T),imgaussfilt(C_norm{ss}(:,T(I)).',[0.0001,2])) 
    caxis([-0.4,0.7])
end


%% hist count of lick events

edges2 = [-2e2:20:3e2];

Lcount = zeros(length(edges2)-1,length(T));
for n = T
    
    M = zeros(length(Pool(n).licktimes),length(Pool(n).licktimes));
    
    for s  = 1:length(Pool(n).licktimes)
        M(s,:) = floor(Pool(n).licktimes(s)*1e3)-floor(Pool(n).licktimes*1e3);
    end
    
    M = reshape(M,[],1);
    
    %     histogram(M,edges)
    countM = histcounts(M,edges2);
    Lcount(:,n) = countM;
    
end

figure
subplot(2,1,2)
shadedErrorBar(edges2(1:end-1),mean(Lcount(:,T),2)/max(max(Lcount(:,T))),...
    std(Lcount(:,T),[],2)/(sqrt(size(Lcount(:,T),2))*max(max(Lcount(:,T)))))
subplot(2,1,1)
for ss = 1:2
    C_norm = zeros(size(C{ss},1),size(C{ss},2));
    for n  = 1:length(C{ss})
        C_norm(:,n) = (C{ss}(:,n)-prctile(C{ss}(:,n),25))/(max(C{ss}(:,n)-prctile(C{ss}(:,n),25))+10);
    end
    
    shadedErrorBar(edges(1:end-1),mean(C_norm(:,T),2),std(C_norm(:,T),[],2)/sqrt(size(C_norm(:,T),2)),'lineProps',c_list{ss})
    hold on
end
      


%% FR of lick-related neurons 

F ={};
window = 10; % ms window
dur = 8000; % length of trial

for ss = 1:4
    F{ss} = zeros(floor(dur/window),length(Pool));
end

P1 = zeros(4,length(Pool));

for n = T
    for ss = 1:4
        for t = 1:floor(dur/window)
            F{ss}(t,n) = mean(mean(rate.PSTH{n,ss}(:,(t-1)*window+1:t*window),2),1);
        end
        % Find peak during stim
        [peak,p_ind] = max(F{ss}(200:300,n));
        spont = mean(F{ss}(1:150,n));
        sd = std(F{ss}(1:150,n),[],1);
        if peak > 4*sd + spont
            if F{ss}(p_ind-1+200,n) > 2*sd + spont && F{ss}(p_ind+1+200,n) > 2*sd + spont
                P1(ss,n) = peak-spont;
            end
        end
        
        
    end
    
    
end
P1 = P1(:,T);


for n = 1:length(T)
    P1(:,n) = P1(:,n)/(max(P1(:,n))+5);
end

sz = 20;
figure
subplot(2,2,1)
scatter(P1(1,:),P1(3,:),sz,'filled')
hold on
yline(0)
xline(0)
axis([-0.1,1.1,-0.1,1.1])

subplot(2,2,2)
scatter(P1(4,:),P1(2,:),sz,'filled')
hold on
yline(0)
xline(0)
axis([-0.1,1.1,-0.1,1.1])

subplot(2,2,3)
scatter(P1(1,:),P1(4,:),sz,'filled')
hold on
yline(0)
xline(0)
axis([-0.1,1.1,-0.1,1.1])

subplot(2,2,4)
scatter(P1(3,:),P1(2,:),sz,'filled')
hold on
yline(0)
xline(0)
axis([-0.1,1.1,-0.1,1.1])

length(unique([find(P1(2,:)>0), find(P1(3,:))]))


    






figure
c_list = {'r','g','b','k'};
for ss = 1:4
    for n = T
        F_norm(:,n) = (F{1,ss}(:,n))/(max(F{1,ss}(:,n))+.5);
    end
    meanF = smooth(mean(F_norm(:,T),2),100);
    sdF = smooth(std(F_norm(:,T),[],2)/sqrt(size(F_norm(:,T),2)),100);
    hold on
%     plot(meanF)
%     pause
    shadedErrorBar(1:8000,meanF,sdF,'lineProps',c_list{ss})
end
        
xlim([1900,3000])
       
%% 06242024 Align to nth lick

edges = [-2e2:5:3e2];

C2 = {};


n= 155; 

for n = T
    event_times= Pool(n).eventtimes(Pool(n).xb.trial_type(:,4)< 5);
    % event_times = event_times + 1.5; % add 1.5 seconds to eventtimes to align to reward period.
    Rt = zeros(length(Pool(n).xb.hit_code),10);
    Spks = {};
    for tr  = 1:length(event_times)
        rt_ind = find(Pool(n).licktimes > event_times(tr) + 1.5 & Pool(n).licktimes < (event_times(tr) + 1.5+2) );
        sp_ind = find(Pool(n).spiketimes > event_times(tr) & Pool(n).spiketimes < (event_times(tr)+2) );
        if length(rt_ind) >10
            Rt(tr,:) = rt_ind(1:10);
        else
            Rt(tr,1:length(rt_ind)) = rt_ind;
        end
        Spks{tr,1} = Pool(n).spiketimes(sp_ind);% -event_times(tr);
    end
    
    
    FR = {};
    for ln = 1:5
        FR{ln}.Hit = [];
    end
    
    for ln = 1:5
        hit_ind = find(Pool(n).xb.hit_code(:,1)==1);
        
        for r = 1:length(hit_ind)
            try
                FR{ln}.Hit = [FR{ln}.Hit;Spks{hit_ind(r)}-Pool(n).licktimes(Rt(hit_ind(r),ln)-1)];
            catch
            end
        end
    end
    %
    % histogram(FR{3}.Hit*1e3,edges)
    % figure
    for ln = 1:5
        countM = histcounts(FR{ln}.Hit*1e3,edges);
        C2{ln}(:,n) = countM;
        %     plot(edges(1:end-1),countM)
        %     hold on
    end
end

% Normalizing C
for n =T
    maxC2 = [];
    for r = 1:5
        maxC2 = [maxC2, max(C2{r}(:,n))];
    end
    maxC2 = max(maxC2);
    for r = 1:5
        C2norm{r}(:,n) = (C2{r}(:,n)-prctile(C2{r}(:,n),25))/(maxC2-prctile(C2{r}(:,n),25)+2);
    end
end

figure
c_list = {'r','g','b','c','k'};

for r = 1:5
    y = mean(C2norm{r}(:,T),2);
    sd = std(C2norm{r}(:,T),[],2)/sqrt(length(T));
    
    
    shadedErrorBar(edges(1:end-1),y,sd,'lineProps',c_list{r})
    hold on
end
    


r = 1;
s_list = zeros(1,length(T));
for t = 1:length(T)
    [~,s_list(t)] = max(C2norm{r}(:,T(t)));
end
[~,I] = sort(s_list);

figure
imagesc(edges(1:end-1),1:length(T),imgaussfilt(C2norm{r}(:,T(I)).',[0.0001,2])) 
% caxis([-0.4,0.7])


























