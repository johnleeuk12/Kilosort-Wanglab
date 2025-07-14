
for n = 1:length(Pool)
%     [~,arg] = min(abs(dataset{n,3}(:,2)-dataset{n,4}(1,2)));
    arg = dataset{n,5} + 200;
    Pool(n).ttr = arg;
end



%{

Edit log 
=======================================

2025/02/10 JHL

Calculating gain modulation for manuscript figure 1
Currently the code isn't optimized to be run as is
after running "determining gain" for either IC or AC data, 
rename gain and stim variables into gain_IC and stim_IC or 
gain_AC and stim_AC respectively before running the rest of the code



%}

% % 
% n_arr = 1:length(Pool);
% c_list = {'r-','b--','r--','b-'};
% % l_list = {'solid','dotted','dotted','solid'};
% list_n = n_arr(list3);
% x_axis = -2:1e-3:6;
% x_axis = x_axis(1:end-1);
% for n = 1:length(list_n)
%     nn = list_n(n);
%     figure(n)
%     for pp = 1:4
%         y = smoothdata(mean(rate.PSTH{nn,pp},1),"gaussian",50);
%         plot(x_axis,y,c_list{pp},'LineWidth',2)
%         hold on
%     end
%     xlim([-0.5,1.5])
%     drawnow
%     pause(0.1)
% end


%% determining gain
clear raster2 rate Lick2 Lick_rate
addpath('D:\GitHub\Kilosort-Wanglab\Analysis_postphy_core')
% modify Pool

for n  = 1:length(Pool)
    for tr  = 1: length(Pool(n).xb.trial_type)
        if tr <201
            Pool(n).xb.trial_type(tr,5) = Pool(n).xb.trial_type(tr,4);
        elseif tr <=460 && tr > Pool(n).ttr
            Pool(n).xb.trial_type(tr,5) = Pool(n).xb.trial_type(tr,4);
        elseif tr >200 && tr <= Pool(n).ttr
            if mod(Pool(n).xb.trial_type(tr,4),2) == 1
                Pool(n).xb.trial_type(tr,5) = 7;
                
            else
                Pool(n).xb.trial_type(tr,5) = 8;
            end
            
        else
            Pool(n).xb.trial_type(tr,5) = Pool(n).xb.trial_type(tr,4);
        end
        if tr >Pool(n).ttr && tr <= 260
            if mod(Pool(n).xb.trial_type(tr,4),2) == 1
                Pool(n).xb.trial_type(tr,5) = 9;
            else
                Pool(n).xb.trial_type(tr,5) = 10;
            end
        end
    end
    for tr_type = 1:10
        [~,ind] = sort(find(Pool(n).xb.trial_type(:,5)==tr_type));
        Pool(n).xb.trial_type(find(Pool(n).xb.trial_type(:,5)==tr_type),6) = ind;
    end
    
end
       
            

% [raster2, rate, Lick2, Lick_rate] = gather_raster_ephys(2, 6, Pool);

[raster2, rate, Lick2, Lick_rate] = gather_raster_ephys2(2, 6, Pool);


%%

list_n = 1:length(Pool);
% list_n = [1:102,126:179];
% list_n = list_n([Pool(:).ttr] > 206);
gain = {};
stim = {};
supp = {};
for r = 1:2
    gain{r}.onset = zeros(1,length(list_n));
    gain{r}.sus = zeros(1,length(list_n));
    gain{r}.all =  zeros(1,length(list_n));
    stim{r}.onset = zeros(1,length(list_n));
    stim{r}.sus = zeros(1,length(list_n));
    stim{r}.all = zeros(1,length(list_n));
    gain{r}.pv =  zeros(1,length(list_n));
    supp{r}.sus =  zeros(1,length(list_n));
end
alpha = 1;
for n = 1:length(list_n)
    nn = list_n(n);
    for r = 1:2 % 5 or 10 kHz, analysed separately
        if r == 1
            y_go =  smoothdata(mean(rate.PSTH{nn,1},1),"gaussian",50);
            raw_go = rate.PSTH{nn,1};
            y_ng = smoothdata(mean(rate.PSTH{nn,3},1),"gaussian",50);
            raw_ng = rate.PSTH{nn,3};
        else
            y_go =  smoothdata(mean(rate.PSTH{nn,4},1),"gaussian",50);
            raw_go = rate.PSTH{nn,4};
            y_ng = smoothdata(mean(rate.PSTH{nn,2},1),"gaussian",50);
            raw_ng = rate.PSTH{nn,2};
            
        end
        
        
        if (mean(y_go(2000:2100)) > mean(y_go(500:1500))+2*std(mean(raw_go(:,500:1500),2))) ...
                || (mean(y_ng(2000:2100)) > mean(y_ng(500:1500))+2*std(mean(raw_ng(:,500:1500),2)))
            stim{r}.onset(n) = 1;
            gain{r}.onset(n) = ((mean(y_go(2000:2100))- mean(y_go(500:1500)))- ...
                (mean(y_ng(2000:2100))- mean(y_ng(500:1500))))/ ...
                (max([(mean(y_go(2000:2100))- mean(y_go(500:1500))), ...
                (mean(y_ng(2000:2100))- mean(y_ng(500:1500)))]) + ...
                alpha);
        end
        
        if (mean(y_go(2100:2450)) > mean(y_go(500:1500))+2*std(mean(raw_go(:,500:1500),2))) ...
                || (mean(y_ng(2100:2450)) > mean(y_ng(500:1500))+2*std(mean(raw_ng(:,500:1500),2)))
            stim{r}.sus(n) = 1;
            gain{r}.sus(n) = ((mean(y_go(2000:2450))- mean(y_go(500:1500)))- ...
                (mean(y_ng(2100:2450))- mean(y_ng(500:1500))))/ ...
                (max([(mean(y_go(2100:2450))- mean(y_go(500:1500))), ...
                (mean(y_ng(2100:2450))- mean(y_ng(500:1500)))]) + ...
                alpha);
        end
        
        % suppression during stim
        if (mean(y_go(2100:2450)) < mean(y_go(500:1500)))
            d1 = y_go(2100:2450);
            d2  = y_go(500:1500);
            bootstat = bootstrp(100,@mean,d1);
            bootstat2 = bootstrp(100,@mean,d2);
            [~,p] = kstest2(bootstat,bootstat2);
            if p < 0.01
                if stim{r}.onset(n) == 0
                    if mean(y_go(2100:2450)) < mean(y_go(500:1500)) - 2*std(y_go(500:1500))
                        supp{r}.sus(n) = 1;
                    end
                end
            end
        elseif (mean(y_ng(2100:2450)) < mean(y_ng(500:1500)))
            d1 = y_ng(2100:2450);
            d2  = y_ng(500:1500);
            bootstat = bootstrp(100,@mean,d1);
            bootstat2 = bootstrp(100,@mean,d2);
            [~,p] = kstest2(bootstat,bootstat2);
            if p < 0.01
                if stim{r}.onset(n) == 0
                    if mean(y_ng(2100:2450)) < mean(y_ng(500:1500)) - 2*std(y_ng(500:1500))
                        supp{r}.sus(n) = 1;
                    end
                end
            end
        end
        
        % remove inhibitory responses
        if max([(mean(y_go(2100:2450))- mean(y_go(500:1500))), ...
                (mean(y_ng(2100:2450))- mean(y_ng(500:1500)))])...
                < 0
            gain{r}.sus(n) = 0;
        end
        if max([(mean(y_go(2000:2100))- mean(y_go(500:1500))), ...
                (mean(y_ng(2000:2100))- mean(y_ng(500:1500)))])...
                < 0
            gain{r}.onset(n) = 0;
        end
        % remove units with 0 spikes in most trials during this period
        if median(y_go(2000:2100)) == 0 && median(y_ng(2000:2100)) == 0
            gain{r}.onset(n) = 0;
            stim{r}.onset(n) = 0;
        end
        
        if median(y_go(2100:2450)) == 0 && median(y_ng(2100:2450)) == 0
            gain{r}.sus(n) = 0;
            stim{r}.sus(n) = 0;
        end
        % removing units end
        
        if abs(gain{r}.onset(n)) > abs(gain{r}.sus(n))
            gain{r}.all(n) = gain{r}.onset(n);
            d = (mean(raw_go(:,2000:2100),2)-mean(y_go(500:1500)));
            d2 = (mean(raw_ng(:,2000:2100),2)-mean(y_ng(500:1500)));
            d  = d + 2*max(mean(y_go(500:1500)),mean(y_ng(500:1500)));
            d2 = d2 + 2*max(mean(y_go(500:1500)),mean(y_ng(500:1500)));
            
            bootstat = bootstrp(100,@mean,d);
            bootstat2 = bootstrp(100,@mean,d2);
            %             [~,p] = kstest2(bootstat,bootstat2);
            %             figure
            %             histfit(bootstat,100,'normal')
            %             hold on
            %             histfit(bootstat2,100,'normal')
            %             p = spk_ttest2(bootstat,bootstat2);
            p = spk_ttest2(d,d2);
            %             p = spk_ttest2(d2,d);
            %             [~,p] = ttest2(d,d2);
            %             [~,p] = kstest2((mean(raw_go(:,2000:2100),2)-mean(y_go(500:1500)))...
            %                 ,(mean(raw_ng(:,2000:2100),2)-mean(y_ng(500:1500))));
            gain{r}.pv(n) = p;
        elseif abs(gain{r}.onset(n)) < abs(gain{r}.sus(n))
            gain{r}.all(n) = gain{r}.sus(n);
            d = (mean(raw_go(:,2100:2450),2)-mean(y_go(500:1500)));
            d2 = (mean(raw_ng(:,2100:2450),2)-mean(y_ng(500:1500)));
            d  = d + 2*max(mean(y_go(500:1500)),mean(y_ng(500:1500)));
            d2 = d2 + 2*max(mean(y_go(500:1500)),mean(y_ng(500:1500)));
            bootstat = bootstrp(100,@mean,d);
            bootstat2 = bootstrp(100,@mean,d2);
            %             [~,p] = kstest2(bootstat,bootstat2);
            p = spk_ttest2(d,d2);
            %             [~,p] = kstest2((mean(raw_go(:,2100:2450),2)-mean(y_go(500:1500)))...
            %                 ,(mean(raw_ng(:,2100:2450),2)-mean(y_ng(500:1500))));
            gain{r}.pv(n) = p;
        end
        
    end
end

for r = 1:2
    stim{r}.all = ((stim{r}.onset + stim{r}.sus)>0);
end

supp_all = (supp{1}.sus | supp{2}.sus);
sum(supp_all& not(stim{1}.all|stim{2}.all))

%% 
addpath('D:\GitHub\Kilosort-Wanglab\Analysis_postphy\Violinplot-Matlab-master');

for r = 1:2
    p_list{r} = (0 < gain{r}.pv) & (gain{r}.pv < 0.05);
   
end

for r = 1:2
    g_list{r} = (0.05 < gain{r}.pv) | (gain{r}.pv ==0);
    g_list{r} = logical(g_list{r}.*stim{r}.all);

end

vd = [gain{1}.all( p_list{1}),gain{2}.all( p_list{2})];
vc = [ones(sum(p_list{1}),1)*2;ones(sum(p_list{2}),1)*1];
vc = vc.';

vd2 = [gain{1}.all( g_list{1}),gain{2}.all( g_list{2})];
vc2 = [ones(sum(g_list{1}),1)*2;ones(sum(g_list{2}),1)*1];
vc2 = vc2.';


color = {[0,0,1],[1,0,1],[0,0,1],[1,0,1]};
figure
% vs2 = violinplot(vd2,vc2,'Orientation', 'horizontal','ViolinColor',[0,0,0]);
hold on
vs1 = violinplot(vd,vc,'Orientation', 'horizontal','ViolinColor' ,[[1,0,1];[0,0,1]] );
hold off
% vs1(1,1).ScatterPlot.MarkerFaceColor = [1,1,1];

xlim([-1.2,1.2])


sum((stim{1}.all+stim{2}.all)>0)
r = 1
median(gain{r}.all(p_list{r}))
std(gain{r}.all(p_list{r}))/sqrt(sum(p_list{r}))

[~,p ] = kstest(gain{r}.all(p_list{r}))



%% gather gain

gain_all = {};
gain_all{1}{1} = gain;
gain_all{1}{2} = p_list;

%% comparing RT and R2 
% sum(gain_all{1}{2}{1}.*gain_all{1}{2}{2})

gain_total = zeros(4,length(Pool));
gain_total(1,:) = gain_all{1}{1}{1}.all;
gain_total(2,:) = gain_all{1}{1}{2}.all;
gain_total(3,:) = gain_all{2}{1}{1}.all;
gain_total(4,:) = gain_all{2}{1}{2}.all;

p_list_total = zeros(8,length(Pool));
p_list_total(1,:) = gain_all{1}{2}{1};
p_list_total(2,:) = gain_all{1}{2}{2};
p_list_total(3,:) = gain_all{2}{2}{1};
p_list_total(4,:) = gain_all{2}{2}{2};
p_list_total(5,:) = (p_list_total(1,:) | p_list_total(3,:));
p_list_total(6,:) = (p_list_total(2,:) | p_list_total(4,:));
p_list_total(7,:) = (p_list_total(1,:) & not(p_list_total(3,:)));
p_list_total(8,:) = (p_list_total(2,:) & not(p_list_total(4,:)));
p_list_total = logical(p_list_total);


%% Code comparing AC and IC
vd = [gain_total(2,p_list_total(4,:)),gain_total(4,p_list_total(4,:))];
vc = [ones(sum(p_list_total(4,:)),1)*1;ones(sum(p_list_total(4,:)),1)*2];
vc = vc.';
figure
vs1 = violinplot(vd,vc,'Orientation', 'horizontal','ViolinColor' ,[1,0,1], ...
                'ShowData', true,'Orientation','vertical' );
% hold on
% scatter(vc,vd)
% for sc = 1:length(vd)/2
%     plot([1,2],[vd(sc),vd(sc+length(vd)/2)])
% end
ylim([-1,1.0])
figure
boxplot(vd,vc)
hold on 
scatter(vc,vd)
for sc = 1:length(vd)/2
    plot([1,2],[vd(sc),vd(sc+length(vd)/2)])
end
ylim([-1,1.0])

grp1 = gain_total(2,p_list_total(4,:));
grp2 = gain_total(4,p_list_total(4,:));

signrank(grp1,grp2)
[~,p] = ttest(grp1,grp2)

median(grp1)
std(grp1)/sqrt(length(grp1))

[~,p ] = kstest(grp1)

figure
vd = [gain_total(2,p_list_total(8,:)),gain_total(4,p_list_total(8,:))];
vc = [ones(sum(p_list_total(8,:)),1)*1;ones(sum(p_list_total(8,:)),1)*2];
vc = vc.';
% vs1 = violinplot(vd,vc,'Orientation', 'horizontal','ViolinColor' ,[1,0,1],'ShowData', false );
boxplot(vd,vc)
hold on
scatter(vc,vd)
ylim([-1,1.5])

mean(gain_total(4,p_list_total(8,:)))
std(gain_total(4,p_list_total(8,:)))/length(gain_total(4,p_list_total(8,:)))

[~,p ] = kstest(gain_total(4,p_list_total(8,:)))

%%
% 
% for r = 1:2
%     p_list{r} = (0 < gain_IC2{r}.pv) & (gain_IC2{r}.pv < 0.01);
%     p_list2{r} = (0 < gain_AC2{r}.pv) & (gain_AC2{r}.pv < 0.01);
%    
% end
% 
% edges= -1.5:0.1:1.5;
% for r = 1:2
%     figure(r+2)
%     histogram(gain_IC2{r}.all( p_list{r}),edges)
% %     hold on
% %     histogram(gain_AC2{r}.all( p_list2{r}),edges)
%     ylim([0,20])
%     hold off
% end
% 
% for r = 1:2
%     figure(r)
%     x = [gain_AC{r}.all( p_list2{r}),gain_IC{r}.all( p_list{r})];
%     x = x.';
%     g = [ones(sum(p_list2{r}),1);ones(sum(p_list{r}),1)*2];
%     boxplot(x,g)
%     hold on 
%     scatter(g,x,20,'k','filled')
%     ylim([-1,1.5])
%     hold off
% end
%     

%%

for r = 1:2
    p_list{r} = (0 < gain_IC{r}.pv) & (gain_IC{r}.pv < 0.05);
    p_list2{r} = (0 < gain_AC{r}.pv) & (gain_AC{r}.pv < 0.05);
   
end

for r = 1:2
    g_list{r} = (0.05 < gain_IC{r}.pv) | (gain_IC{r}.pv ==0);
    g_list{r} = logical(g_list{r}.*stim_IC{r}.all);
    g_list2{r} = (0.05 < gain_AC{r}.pv) | (gain_AC{r}.pv ==0);
    g_list2{r} = logical(g_list2{r}.*stim_AC{r}.all);
end



edges= -1.5:0.1:1.5;
for r = 1:2
    figure(r)
%     histogram(gain_AC{r}.all( p_list2{r}),edges)
    hold on
    histogram(gain_IC{r}.all( p_list{r}),edges)
    xline(0,'linewidth',2)
%     xline(mean(gain_IC{r}.all( p_list{r})),'--r','linewidth',2)
    [f,xi] = ksdensity(gain_IC{r}.all( p_list{r}));
    plot(xi,f*max(histcounts(gain_IC{r}.all( p_list{r}),edges)),'linewidth',2)
    [~,ind] = max(f);
    xline(xi(ind),'--r','linewidth',2)
    ylim([0,15])
    hold off
end

for r = 1:2
    figure(r+2)
    histogram(gain_AC{r}.all( p_list2{r}),edges)
    hold on
%     histogram(gain_IC{r}.all( p_list{r}),edges)
    xline(0,'linewidth',2)
%     xline(mean(gain_AC{r}.all( p_list2{r})),'--r','linewidth',2)
    [f,xi] = ksdensity(gain_AC{r}.all( p_list2{r}));
    plot(xi,f*max(histcounts(gain_AC{r}.all( p_list2{r}),edges)),'linewidth',2)
    [~,ind] = max(f);
    xline(xi(ind),'--r','linewidth',2)
    ylim([0,20])
    hold off
end

%%
% for r = 1:2
%     figure(r)
% %     histogram(gain_AC{r}.all( p_list2{r}),edges)
%     hold on
%     histogram(gain_IC{r}.all( p_list{r}),edges)
%     xline(0,'linewidth',2)
% %     xline(mean(gain_IC{r}.all( p_list{r})),'--r','linewidth',2)
%     [f,xi] = ksdensity(gain_IC{r}.all( p_list{r}));
%     plot(xi,f*max(histcounts(gain_IC{r}.all( p_list{r}),edges)),'linewidth',2)
%     [~,ind] = max(f);
%     xline(xi(ind),'--r','linewidth',2)
%     ylim([0,15])
%     hold off
% end
% 
% for r = 1:2
%     figure(r+2)
%     histogram(gain_AC{r}.all( p_list2{r}),edges)
%     hold on
% %     histogram(gain_IC{r}.all( p_list{r}),edges)
%     xline(0,'linewidth',2)
% %     xline(mean(gain_AC{r}.all( p_list2{r})),'--r','linewidth',2)
%     [f,xi] = ksdensity(gain_AC{r}.all( p_list2{r}));
%     plot(xi,f*max(histcounts(gain_AC{r}.all( p_list2{r}),edges)),'linewidth',2)
%     [~,ind] = max(f);
%     xline(xi(ind),'--r','linewidth',2)
%     ylim([0,20])
%     hold off
% end

vd = [gain_IC{1}.all( p_list{1}),gain_IC{2}.all( p_list{2}),...
    gain_AC{1}.all( p_list2{1}),gain_AC{2}.all( p_list2{2})];
vc = [ones(sum(p_list{1}),1)*4;ones(sum(p_list{2}),1)*3;...
    ones(sum(p_list2{1}),1)*2;ones(sum(p_list2{2}),1)*1];
vc = vc.';

vd2 = [gain_IC{1}.all( g_list{1}),gain_IC{2}.all( g_list{2}),...
    gain_AC{1}.all( g_list2{1}),gain_AC{2}.all( g_list2{2})];
vc2 = [ones(sum(g_list{1}),1)*4;ones(sum(g_list{2}),1)*3;...
    ones(sum(g_list2{1}),1)*2;ones(sum(g_list2{2}),1)*1];
vc2 = vc2.';


color = {[0,0,1],[1,0,1],[0,0,1],[1,0,1]};
figure
vs2 = violinplot(vd2,vc2,'Orientation', 'horizontal','ViolinColor',[0,0,0]);
hold on
vs1 = violinplot(vd,vc,'Orientation', 'horizontal','ViolinColor' ,[[1,0,1];[0,0,1];[1,0,1];[0,0,1]] );
hold off
% vs1(1,1).ScatterPlot.MarkerFaceColor = [1,1,1];

xlim([-1.2,1.2])

% 
% for r = 1:2
%     figure(r)
%     x = [gain_AC{r}.all( p_list2{r}),gain_IC{r}.all( p_list{r})];
%     x = x.';
%     g = [ones(sum(p_list2{r}),1);ones(sum(p_list{r}),1)*2];
%     boxplot(x,g)
%     hold on 
%     scatter(g,x,20,'k','filled')
%     ylim([-1,1.5])
%     hold off
% end
%     
% 
% 
% r = 2
% [~,p ] = kstest2(gain_IC{r}.all(p_list{r}),gain_AC{r}.all(p_list2{r}))


%% 
n_list = 1:length(Pool);
for n = n_list(p_list2{2})
    figure(n)
    for pp = [1:4]
        y = smoothdata(mean(rate.PSTH{n,pp},1),"gaussian",50);
        plot(x_axis,y,c_list{pp},'LineWidth',2)
        hold on
    end
    xlim([-0.5,1.5])
    drawnow
    pause(0.1)
end
    

%% 
sum((stim{1}.all+stim{2}.all)>0)
r = 1
median(gain_IC{r}.all(p_list{r}))
std(gain_IC{r}.all(p_list{r}))/sqrt(sum(p_list{r}))

[~,p ] = kstest(gain_AC{r}.all(p_list2{r}))



