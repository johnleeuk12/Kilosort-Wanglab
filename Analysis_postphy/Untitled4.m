 
% 
n_arr = 1:length(Pool);
c_list = {'r-','b--','r--','b-'};
% l_list = {'solid','dotted','dotted','solid'};
list_n = n_arr(list3);
x_axis = -2:1e-3:6;
x_axis = x_axis(1:end-1);
for n = 1:length(list_n)
    nn = list_n(n);
    figure(n)
    for pp = 1:4
        y = smoothdata(mean(rate.PSTH{nn,pp},1),"gaussian",50);
        plot(x_axis,y,c_list{pp},'LineWidth',2)
        hold on
    end
    xlim([-0.5,1.5])
    drawnow
    pause(0.1)
end


%% determining gain
[raster2, rate, Lick2, Lick_rate] = gather_raster_ephys(2, 6, Pool);

list_n = 1:length(Pool);
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
alpha = 2;
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
                if stim{r}.onset(n) == 1
                    supp{r}.sus(n) = 1;
                end
            end
        elseif (mean(y_ng(2100:2450)) < mean(y_ng(500:1500)))
            d1 = y_ng(2100:2450);
            d2  = y_ng(500:1500);
            bootstat = bootstrp(100,@mean,d1);
            bootstat2 = bootstrp(100,@mean,d2);
            [~,p] = kstest2(bootstat,bootstat2);
            if p < 0.01
                if stim{r}.onset(n) == 1
                    supp{r}.sus(n) = 1;
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
    p_list{r} = (0 < gain_IC{r}.pv) & (gain_IC{r}.pv < 0.10);
    p_list2{r} = (0 < gain_AC{r}.pv) & (gain_AC{r}.pv < 0.10);
   
end

for r = 1:2
    g_list{r} = (0.10 < gain_IC{r}.pv) | (gain_IC{r}.pv ==0);
    g_list{r} = logical(g_list{r}.*stim_IC{r}.all);
    g_list2{r} = (0.10 < gain_AC{r}.pv) | (gain_AC{r}.pv ==0);
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
r =1     
[~,p ] = kstest(gain_AC{r}.all(p_list{r}))
mean(gain_AC{r}.all(p_list{r}))



