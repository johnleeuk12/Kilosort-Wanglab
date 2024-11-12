% Calculating times for lick_bout onsets and offsets


n = 155;
ILI = Pool(n).licktimes(2:end)-Pool(n).licktimes(1:end-1);
lick_onset= find(ILI>2); % lick bout = 2 seconds
lick_onset = lick_onset+1;
lick_offset = lick_onset-1;

event_times= Pool(n).eventtimes(Pool(n).xb.trial_type(:,4)< 5);



L_all = zeros(1,floor((Pool(n).licktimes(end)+8)*1e3)); % in ms, entire array of licks converted to binary
L_all(floor(Pool(n).licktimes(lick_offset)*1e3)) = 1;

L_mat = zeros(length(event_times),8*1e3);




pre = 1;
post = 7;

for e = 1:length(event_times)
    L_mat(e,:) = L_all(floor((event_times(e)-pre)*1e3):floor((event_times(e)+post)*1e3)-1);
end


plot(smooth(mean(L_mat,1)))



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





%% 06242024 Align to nth lick

bin = 10;
edges = -1e3:bin:1e3;


pre = 2;
post = 6;

% n= 155; 
w = gausswin(5);

C2 = {};
C2.hit = zeros(length(Pool),length(edges)-1);
C2.hit =C2.hit.';
C2.FA = zeros(length(Pool),length(edges)-1);
C2.FA =C2.FA.';

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
        event_times= Pool(n).eventtimes(Pool(n).xb.trial_type(1:length(Pool(n).eventtimes),4)< 5);
    catch
        event_times= Pool(n).eventtimes(Pool(n).xb.trial_type(:,4)< 5);
    end
    ILI = Pool(n).licktimes(2:end)-Pool(n).licktimes(1:end-1);
    lick_onset= find(ILI>2); % lick bout = 2 seconds
    lick_onset = lick_onset+1;
    lick_offset = lick_onset-1;
    
    % event_times = event_times + 1.5; % add 1.5 seconds to eventtimes to align to reward period.
    Rt = zeros(length(Pool(n).xb.hit_code),10);
    Rt2 = zeros(length(Pool(n).xb.hit_code),10);
    spont_raw = zeros(length(event_times),1);
    Spks = {};
    Lcks = {};
    for tr  = 1:length(event_times)-1
        lick_ind =  find(Pool(n).licktimes > event_times(tr) & ...
                        Pool(n).licktimes < (event_times(tr+1)) ); 
        
        rt_ind = find(Pool(n).licktimes(lick_onset) > event_times(tr) & ...
                        Pool(n).licktimes(lick_onset) < (event_times(tr) + post) );
        rt_ind2 = find(Pool(n).licktimes(lick_onset) > event_times(tr)+3 & ...
                        Pool(n).licktimes(lick_onset) < (event_times(tr+1)) );
                    
        sp_ind = find(Pool(n).spiketimes > event_times(tr) & ...
                        Pool(n).spiketimes < (event_times(tr)+post) );
                    
        sp_ind2 = find(Pool(n).spiketimes > event_times(tr)-2 & ...
                        Pool(n).spiketimes < (event_times(tr)-0.5) ); % spont rate, 1 sec before event
        if length(rt_ind) >10
            Rt(tr,:) = rt_ind(1:10);
        else
            Rt(tr,1:length(rt_ind)) = rt_ind;
        end
        if length(rt_ind2) >10
            Rt2(tr,:) = rt_ind2(1:10);
        else
            Rt2(tr,1:length(rt_ind2)) = rt_ind2;
        end
        Spks{tr,1} = Pool(n).spiketimes(sp_ind);% -event_times(tr);
        spont_raw(tr,1) = length(Pool(n).spiketimes(sp_ind2));
        Lcks{tr,1} = Pool(n).licktimes(lick_ind);
    end
    
    FR = {};
    %     for ln = 1:5
%     ln = 1;
    FR.Hit = [];
    FR.FA = [];
    FR.slick = [];
    LR = {};
    LR.Hit = [];
    LR.FA = [];
    LR.slick = [];
    %     end
    
    %     for ln = 1:5
    hit_ind = find(Pool(n).xb.hit_code(:,1)==1);
    FA_ind = find(Pool(n).xb.hit_code(:,3)==1);

  
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
    for r =1:length(Rt2)
        if Rt2(r,1) >0
            LR.slick = [LR.slick; Lcks{r,1}-Pool(n).licktimes(lick_onset(Rt2(r,1)))];
        end
        try
            FR.slick = [FR.slick;Spks{r}-Pool(n).licktimes(lick_onset(Rt2(r,1)))];
        catch
        end
        try
            FR.slick = [FR.slick;Spks{r}-Pool(n).licktimes(lick_onset(Rt2(r,2)))];
        catch
        end
    end
    
    
    nb_hit2 = sum(sum(Rt(hit_ind,:)>0)); 
    nb_FA2 = sum(sum(Rt(FA_ind,:)>0));
     nb_slick = sum(sum(Rt2(:,:)>0));
    %     end
    
    %
%     histogram(FR{1}.Hit*1e3,edges)
    %     % figure
    
    countM = histcounts(FR.Hit*1e3,edges);
%     C2.hit(:,n) = countM*1e3/(bin*size(hit_ind,1)); % 5 ms bin, converting to FR 
    C2.hit(:,n) = countM*1e3/(bin*nb_hit2);
    countM = histcounts(FR.FA*1e3,edges);
%     C2.FA(:,n) = countM*1e3/(bin*size(FA_ind,1));
    C2.FA(:,n) = countM*1e3/(bin*nb_FA2);
        countM = histcounts(FR.slick*1e3,edges);
%     C2.FA(:,n) = countM*1e3/(bin*size(FA_ind,1));
    C2.slick(:,n) = countM*1e3/(bin*nb_slick);
    spont.mean(n,1) = mean(spont_raw); % in FR, since we're counting nb of spks per sec
    spont.mean(n,2) = mean(spont_raw(hit_ind));
%     spont.mean(n,3) = mean(spont_raw(FA_ind));
    spont.std(n,1) = std(spont_raw,[],1);
    spont.std(n,2) = std(spont_raw(hit_ind),[],1);
%     spont.std(n,3) = std(spont_raw(FA_ind),[],1);
%     figure(n)
%     plot(edges(1:end-1),filter(w,1,countM))
%     drawnow()
%     pause(0.1)
    L.hit(:,n) = histcounts(LR.Hit*1e3,edges)*1e3/(bin*nb_hit2);
    L.FA(:,n) = histcounts(LR.FA*1e3,edges)*1e3/(bin*nb_FA2);
    L.slick(:,n) = histcounts(LR.slick*1e3,edges)*1e3/(bin*nb_slick);
%     
    
end

%%

% Normalizing C
C2norm = {};
C2norm.hit = zeros(length(edges)-1,length(Pool)); 
C2norm.FA = zeros(length(edges)-1,length(Pool)); 
C2norm.slick = zeros(length(edges)-1,length(Pool)); 

% C2norm = zeros(length(edges),length(Pool));

w = gausswin(5); %gaussian smoothing
C3 = {};
% C3.hit = C2.hit;
% C3.FA = C2.FA;
% C3.hit = imfilter(C2.hit,w,'replicate');
% C3.FA = imfilter(C2.FA,w,'replicate');
C3.hit = smoothdata(C2.hit,"gaussian",20); 
C3.FA = smoothdata(C2.FA,"gaussian",20);
C3.slick = smoothdata(C2.slick,"gaussian",20);

% C3.hit = imgaussfilt(C2.hit,[2,1]);
% C3.FA = imgaussfilt(C2.FA,[2,1]);
% C3.hit = C2.hit;
% C3.FA = C2.FA;


max_ind = {};
max_ind.hit = [];
max_ind.FA = [];
max_ind.slick = [];
thresh = zeros(1,length(Pool));
% m_ind = {};
for n =1:length(Pool)
    [~,m1] = max(C3.hit(:,n));
    [~,m2] = max(C3.FA(:,n));
    [~,m3] = max(C3.slick(:,n));
%     try
%         maxC2 = max([mean(C3.hit(m1-1:m1+1,n)),mean(C3.hit(m1-1:m1+1,n))]);
%     catch
        maxC2 = max([mean(C3.hit(m1,n)),mean(C3.FA(m2,n))]);
%     end
    %     maxC2 = max([mean(C3.hit(m1-1:m1+1,n)),mean(C3.hit(m1-1:m1+1,n))]);
    %     C2norm(:,n) = (C3(:,n)-prctile(C3(:,n),25))/(maxC2-prctile(C3(:,n),25)+2);
    C2norm.hit(:,n) = (C3.hit(:,n)-spont.mean(n,1))/(abs((maxC2-spont.mean(n,1)))+.5);
    %     C2norm.hit(:,n) = (C3.hit(:,n)-spont.hit(n,1))/(maxC2+5);
    
    C2norm.FA(:,n) = (C3.FA(:,n)-spont.mean(n,1))/(abs((maxC2-spont.mean(n,1)))+.5);
    %     C2norm.FA(:,n) = (C3.FA(:,n)-spont.FA(n,1))/(maxC2+5);
%     thresh(1,n) = spont.std(n,1)/(maxC2-spont.mean(n,1)+5);
    C2norm.slick(:,n) = (C3.slick(:,n)-spont.mean(n,1))/(abs((maxC2-spont.mean(n,1)))+.5);

    try
        if mean(C3.hit(m1-1:m1+1,n))> 2*spont.std(n,1)+spont.mean(n,1)
            max_ind.hit = [max_ind.hit, n];
        end
        if max(C3.FA(m2-1:m2+1,n))> 2*spont.std(n,1)+spont.mean(n,1)
            max_ind.FA = [max_ind.FA, n];
        end
        if max(C3.slick(m2-1:m2+1,n))> 2*spont.std(n,1)+spont.mean(n,1)
            max_ind.slick = [max_ind.slick, n];
        end
    catch
    end
end




max_ind.all = intersect(max_ind.hit,max_ind.FA);


figure
c_list = {'r','g','b','c','k'};
% C2norm2 = C2norm(:,max_ind);
% C2norm.hit = C3.hit-spont.mean.' ;
% C2norm.FA = C3.FA-spont.mean.' ;


y = mean(C2norm.hit(:,max_ind.hit),2);

ys = [];
x = edges(1:end-1).';
for i = 1:length(y)
    ys(i) = gaussian_kern_reg(x(i),x,y,5);
end

ys = smoothdata(y,"gaussian",20);

sd = std(C2norm.hit(:,max_ind.hit),[],2)/sqrt(size(max_ind.hit,2));
% shadedErrorBar(edges(1:end-1),movmean(y,5,"Endpoints","shrink"),sd,'lineProps',c_list{3})
shadedErrorBar(edges(1:end-1),y,sd,'lineProps',c_list{3})
hold on
y = mean(C2norm.FA(:,max_ind.FA),2);
ys = [];
x = edges(1:end-1).';
for i = 1:length(y)
    ys(i) = gaussian_kern_reg(x(i),x,y,5);
end
ys = smoothdata(y,"gaussian",20);
sd = std(C2norm.FA(:,max_ind.FA),[],2)/sqrt(size(max_ind.FA,2));
% shadedErrorBar(edges(1:end-1),movmean(y,5,"Endpoints","shrink"),sd,'lineProps',c_list{1})
shadedErrorBar(edges(1:end-1),y,sd,'lineProps',c_list{1})

y = mean(C2norm.slick(:,max_ind.slick),2);
x = edges(1:end-1).';
sd = std(C2norm.slick(:,max_ind.slick),[],2)/sqrt(size(max_ind.slick,2));
% shadedErrorBar(edges(1:end-1),movmean(y,5,"Endpoints","shrink"),sd,'lineProps',c_list{1})
shadedErrorBar(edges(1:end-1),y,sd,'lineProps','g')




cmap = colormap(gray);
set(gca,'ColorScale','log')
caxis([0.0001 0.1])  
for t = 5:length(y)-5
%     [h,p] = ttest2(reshape(C2norm.hit(t-2:t+2,max_ind.hit),[],1),...
%         reshape(C2norm.FA(t-2:t+2,max_ind.FA),[],1));
    [~,p] = ttest2(C2norm.hit(t,max_ind.hit),C2norm.FA(t,max_ind.FA));
    if p < 0.05
        scatter(edges(t),0.5,50,p,'|','linewidth',3)
    end
end
  

axis([-1000,1000,-0.1,0.6])



figure
imagesc(smoothdata(C2norm.slick(:,max_ind.hit),1).')
caxis([-0.05 1])  
% 
figure
plot(smoothdata(mean(L.slick(:,:),2).',"gaussian",5))

figure
imagesc(smoothdata(C2norm.hit(:,max_ind.hit),1).')
caxis([-0.05 1])  

%%
Av = {};
max_ind.all2 = unique([max_ind.hit, max_ind.FA]);
Av.pre= zeros(length(Pool),2);
Av.post = zeros(length(Pool),2);
win = floor(100/bin);
onset = floor((length(edges)-1)/2);


for n = 1:length(max_ind.all2)
    Av.pre(max_ind.all2(n),1) = mean(C2norm.hit(onset-win:onset,max_ind.all2(n)));
    Av.pre(max_ind.all2(n),2) = mean(C2norm.FA(onset-win:onset,max_ind.all2(n)));
    Av.post(max_ind.all2(n),1) = mean(C2norm.hit(onset+1:onset+win,max_ind.all2(n)));
    Av.post(max_ind.all2(n),2) = mean(C2norm.FA(onset+1:onset+win,max_ind.all2(n)));    
end

list_pre2 = zeros(4,length(Pool));
for n= 1:length(Pool)
    if Av.pre(n,1) > thresh(1,n)
        list_pre2(1,n) = 1;
    end
    if Av.pre(n,2) > thresh(1,n)
        list_pre2(2,n) = 1;
    end
end

list_pre2(3,:) = list_pre2(1,:).*list_pre2(2,:);
list_pre2(4,:) =  (list_pre2(1,:)>0) + (list_pre2(2,:)>0);
list_pre2 = list_pre2.' ;
list_pre2 = logical(list_pre2);


sum(list_pre2)
size(max_ind.all2,2)

figure
[c,S] = polyfit(Av.pre(:,1),Av.pre(:,2),1);
[y_fit,delta] = polyval(c,Av.pre(:,1),S);
R_squared = 1 - (S.normr/norm(Av.pre(:,2) - mean(Av.pre(:,2))))^2;

scatter(Av.pre(:,1),Av.pre(:,2));
hold on
xline([0,0.1,-0.1])
yline([0,0.1,-0.1])
plot(Av.pre(:,1),y_fit,'r--','LineWidth',1)
text(0.5,-0.3,['y = ',num2str(c(1),'%.2f'),'x'])
text(0.5,-0.5,['R^2 = ',num2str(R_squared,'%.2f')])
axis([-1,1,-1,1])

hold off


figure
[c,S] = polyfit(Av.post(:,1),Av.post(:,2),1);
[y_fit,delta] = polyval(c,Av.post(:,1),S);
R_squared = 1 - (S.normr/norm(Av.post(:,2) - mean(Av.post(:,2))))^2;

hold on
scatter(Av.post(:,1),Av.post(:,2));
plot(Av.post(:,1),y_fit,'r--','LineWidth',1)
xline([0,0.1,-0.1])
yline([0,0.1,-0.1])
text(0.5,-0.3,['y = ',num2str(c(1),'%.2f'),'x'])
text(0.5,-0.5,['R^2 = ',num2str(R_squared,'%.2f')])
axis([-1,1,-1,1])

hold off


%%

% Calculate number of neurons with less than 1spk per sec

% bin into 20ms
binsize = 20;
bin2 = floor(binsize/bin);
edges2  =1:bin2:length(edges)-1;
C4 = {};
C4.hit = zeros(length(1:bin2:length(edges)-1),length(Pool));
C4.FA = zeros(length(1:bin2:length(edges)-1),length(Pool));

disp(num2str(sum(mean(C4.hit,1)>1)))


t2 =1;
for t = 1:bin2:length(edges)-1
    C4.hit(t2,:) = mean(C2.hit(t:t+bin2-1,:),1);
    C4.FA(t2,:) = mean(C2.FA(t:t+bin2-1,:),1);
    t2 = t2+1;
end


% 20ms window
win = floor(100/binsize);
onset = floor((size(C4.hit,1)/2));

list_pre = zeros(length(Pool),2);
list_post = zeros(length(Pool),2);
for n = 1:length(Pool)
    if mean(C4.hit(:,n),1)>1 && mean(C4.FA(:,n),1)>0
        if max(C4.hit(onset-win:onset,n)) > spont.mean(n,2) + 2*spont.std(n,2)
            list_pre(n,1) = 1;
        end
        if max(C4.hit(onset+1:onset+win,n)) > spont.mean(n,2) + 2*spont.std(n,2)
            list_post(n,1) = 1;
        end
        
        if max(C4.FA(onset-win:onset,n)) > spont.mean(n,2) + 2*spont.std(n,2)
            list_pre(n,2) = 1;
        end
        if max(C4.FA(onset+1:onset+win,n)) > spont.mean(n,2) + 2*spont.std(n,2)
            list_post(n,2) = 1;
        end
    end
end

list_pre(:,3) = list_pre(:,1).*list_pre(:,2);
list_post(:,3) = list_post(:,1).*list_post(:,2);
list_pre = logical(list_pre);
list_post = logical(list_post);

sum(list_pre)
sum(list_post)

figure
c_list = {'r','g','b','c','k'};

list_pre(:,1) = list_pre2(:,3);

y = mean(C2.hit(:,list_pre(:,1)),2);
ys = [];
x = edges(1:end-1).';
for i = 1:length(y)
    ys(i) = gaussian_kern_reg(x(i),x,y,20);
end
sd = std(C2.hit(:,list_pre(:,1)),[],2)/sqrt(sum(list_pre(:,1)));
shadedErrorBar(x,ys,sd,'lineProps',c_list{3})
hold on
y = mean(C2.FA(:,list_pre(:,1)),2);
ys = [];
x = edges(1:end-1).';
for i = 1:length(y)
    ys(i) = gaussian_kern_reg(x(i),x,y,20);
end
sd = std(C2.FA(:,list_pre(:,1)),[],2)/sqrt(sum(list_pre(:,1)));
shadedErrorBar(edges(1:end-1),ys,sd,'lineProps',c_list{1})

cmap = colormap(gray);
set(gca,'ColorScale','log')
caxis([0.0001 0.1])  
for t = 5:length(y)-5
%     [h,p] = ttest2(reshape(C2norm.hit(t-2:t+2,max_ind.hit),[],1),...
%         reshape(C2norm.FA(t-2:t+2,max_ind.FA),[],1));
    [~,p] = ttest2(C2norm.hit(t,list_pre(:,1)),C2norm.FA(t,list_pre(:,1)));
    if p < 0.05
        scatter(edges(t),0.38,50,p,'|','linewidth',3)
    end
end

figure
sd = std(L.FA,[],2)/sqrt(length(Pool));
y = mean(L.FA,2);
for i = 1:length(y)
    ys(i) = gaussian_kern_reg(x(i),x,y,20);
end
shadedErrorBar(edges(1:end-1),ys,sd,'lineProps',c_list{1})

hold on 
sd = std(L.hit,[],2)/sqrt(length(Pool));
y = mean(L.hit,2);
for i = 1:length(y)
    ys(i) = gaussian_kern_reg(x(i),x,y,20);
end
shadedErrorBar(edges(1:end-1),ys,sd,'lineProps',c_list{3})




% figure
% c_list = {'r','g','b','c','k'};
% % C2norm2 = C2norm(:,max_ind);
% % C2norm.hit = C3.hit-spont.mean.' ;
% % C2norm.FA = C3.FA-spont.mean.' ;
% 
% % C2norm.hit = movmean(C2norm.hit,5,1,'Endpoints','shrink');
% % C2norm.FA = movmean(C2norm.FA,5,1,'Endpoints','shrink');
% y = mean(C4.hit(:,list_post(:,1)),2);
% sd = std(C4.hit(:,list_post(:,1)),[],2)/sqrt(sum(list_post(:,1)));
% shadedErrorBar(edges2,y,sd,'lineProps',c_list{3})
% hold on
% y = mean(C4.FA(:,list_post(:,2)),2);
% sd = std(C4.FA(:,list_post(:,2)),[],2)/sqrt(sum(list_post(:,2)));
% shadedErrorBar(edges2,y,sd,'lineProps',c_list{1})


%% 

post_hit =  mean(C2.hit(300:end,list_pre(:,1)),1);%-mean(C2.hit(25:100,list_pre(:,1)),1);
post_FA =  mean(C2.FA(300:end,list_pre(:,1)),1);%-mean(C2.FA(25:100,list_pre(:,1)),1);

test = post_hit- post_FA;

[T,p] = kstest(test)

figure

bar([1,2],[mean(test),mean(test2)])                

hold on

er = errorbar([1,2],[mean(test),mean(test2)],[std(test)/sqrt(length(test)),std(test2)/sqrt(length(test2))]);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
%%

s_list = {};
% s_list.hit = zeros(1,size(max_ind.all,2));
s_list.FA = zeros(1,sum(list_pre(:,1)));
C5 = C2norm.FA(:,list_pre(:,1));
for n = 1:size(C5,2)
    [~,s_list.FA(1,n)] = max(C5(:,n));
end

[~,I] = sort(s_list.FA);
figure

C6 = C2norm.hit(:,list_pre(:,1));
% subplot(1,2,1)
imagesc(edges2,1:size(C5,2),imgaussfilt(C6(:,I).',[0.0001,5])) 
caxis([-0.0,0.7])
%%






% all units, sorted

% max_ind.all = unique([max_ind.hit,max_ind.FA]);
% 
% C4 = {};
% C4.hit = C2norm.hit(:,max_ind.all);
% C4.FA = C2norm.FA(:,max_ind.all);
% s_list = {};
% s_list.hit = zeros(1,size(max_ind.all,2));
% s_list.FA = zeros(1,size(max_ind.all,2));
% 
% for n = 1:size(max_ind.all,2)
%     nn = max_ind.all(n);
%     [~,s_list.hit(n)] = max(C4.hit(:,n));
%     [~,s_list.FA(n)] = max(C4.FA(:,n));
% end
% 
% I = {};
% [~,I.hit] = sort(s_list.hit);
% [~,I.FA] = sort(s_list.FA);

% figure
% 
% subplot(1,2,1)
% imagesc(edges(1:end-1),1:size(max_ind.all,2),imgaussfilt(C4.hit(:,I.hit).',[0.0001,2])) 
% caxis([-0.3,0.7])
% 
% 
% subplot(1,2,2)
% imagesc(edges(1:end-1),1:size(max_ind.all,2),imgaussfilt(C4.FA(:,I.hit).',[0.0001,2])) 
% caxis([-0.3,0.7])
% title('Hit')
% 
% figure
% subplot(1,2,1)
% imagesc(edges(1:end-1),1:size(max_ind.all,2),imgaussfilt(C4.FA(:,I.FA).',[0.0001,2])) 
% caxis([-0.3,0.7])
% 
% subplot(1,2,2)
% imagesc(edges(1:end-1),1:size(max_ind.all,2),imgaussfilt(C4.hit(:,I.FA).',[0.0001,2])) 
% caxis([-0.3,0.7])
% title('FA')

% %% PCA on these neurons 
% 
% [U,S,V] = svd(C4.FA.');
% 
% K = S*V.';
% 
% 
% figure
% for p = 1:5
%     subplot(1,5,p)
%     plot(K(p,:))
% end
% 
%     
%     
%     
    
    
    
    
    
    
    

