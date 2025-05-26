

% [gain_rt, stim_rt,~,~] = ana_gain(Pool,3)

clearvars -except Pool;

region = 'PPC';
list_r1 = load(['D:\DATA\listTV3_',region,'_R1.mat']);
list_r2 = load(['D:\DATA\listTV3_',region,'_R2.mat']);
list_rt = load(['D:\DATA\listTV3_',region,'_RT.mat']);
list_n = [list_r1.five, list_r2.five, list_rt.five, ...
    list_r1.ten, list_r2.ten, list_rt.ten];


list_np = {};
temp = load(['D:\DATA\listTV3_',region,'_R1p.mat']);
list_np{1,1} = temp.five+1;
list_np{1,2} = temp.ten+1;
temp = load(['D:\DATA\listTV3_',region,'_R2p.mat']);
list_np{2,1} = temp.five+1;
list_np{2,2} = temp.ten+1;
temp = load(['D:\DATA\listTV3_',region,'_RTp.mat']);
list_np{3,1} = temp.five+1;
list_np{3,2} = temp.ten+1;

list_nn = {};
temp = load(['D:\DATA\listTV3_',region,'_R1n.mat']);
list_nn{1,1} = temp.five+1;
list_nn{1,2} = temp.ten+1;
temp = load(['D:\DATA\listTV3_',region,'_R2n.mat']);
list_nn{2,1} = temp.five+1;
list_nn{2,2} = temp.ten+1;
temp = load(['D:\DATA\listTV3_',region,'_RTn.mat']);
list_nn{3,1} = temp.five+1;
list_nn{3,2} = temp.ten+1;




list_n = unique(list_n)+1;
%% 
% [gain_r2, stim_r2,~,~] = ana_gain(Pool,2,1,0.5,list_n);
% [gain_rt, stim_rt,~,~] = ana_gain(Pool,3,1,0.5,list_n);
 
[gain_r2, stim_r2,~,~] = ana_gain3(Pool,2,1,0.5,list_n,list_np,list_nn);
[gain_rt, stim_rt,~,~] = ana_gain3(Pool,3,1,0.5,list_n,list_np,list_nn);
gain_r2{1}.all = -gain_r2{1}.all;
gain_rt{1}.all = -gain_rt{1}.all;
%%
gain_all ={};



% n_list  ={};
p_list_c = {};
for r = 1:2
    p_list_r2{r} = logical(stim_r2{r}.all.*((0 < gain_r2{r}.pv) & (gain_r2{r}.pv < 0.05)));
    p_list_rt{r} = logical(stim_rt{r}.all.*((0 < gain_rt{r}.pv) & (gain_rt{r}.pv < 0.05)));
%     p_list_rt{r} = (0 < gain_rt{r}.pv) & (gain_rt{r}.pv < 0.05);
    p_list_supp{r}= logical(stim_r2{r}.supp.*((0 < gain_r2{r}.pv_supp) & (gain_r2{r}.pv_supp< 0.05)));
    
    p_list_c{r} =((p_list_r2{r}) | (p_list_r2{r}));
    
    gain_all{r} = [gain_r2{r}.all(p_list_c{r});gain_rt{r}.all(p_list_c{r})];
    gain_all{r} = [gain_all{r};zeros(1,length(gain_all{r}))];
    n_list{r} = find(p_list_c{r}==1);
end

% test = [Pool.ttr];
% [~,~,ttr_list]= unique(test);
% 
% n_list
r= 1;
cmap = {[0,0,1],[1,0,1]};
addpath('D:\GitHub\Kilosort-Wanglab\Analysis_postphy\Violinplot-Matlab-master');


for r = 1:2
    figure('Position',[100 100 500 500])
%     edges = -1.3:0.1:1.3;
    edges = -2.0:0.1:2.0;

%     histogram(gain_r2{r}.all(logical(not(p_list_r2{r}).*stim_r2{r}.all)),edges,'facecolor','k')
    subplot(2,1,1)
    histogram(gain_r2{r}.all(stim_r2{r}.all),edges,'facecolor','k')
%     histogram(gain_rt{r}.all(stim_rt{r}.all),edges,'facecolor','k')
    hold on
    histogram(gain_r2{r}.all(p_list_r2{r}),edges,'facecolor',"g")
%     histogram(gain_rt{r}.all(p_list_rt{r}),edges,'facecolor',"g")
    xline(mean(gain_r2{r}.all(stim_r2{r}.all)),'-r')
    xline(median(gain_r2{r}.all(stim_r2{r}.all)),'--r')
%     xline(mean(gain_rt{r}.all(stim_r2{r}.all)),'-r')
%     xline(median(gain_rt{r}.all(stim_rt{r}.all)),'--r')
    hold off
    ylim([0,10])
    [p,~]= signrank(gain_r2{r}.all(p_list_r2{r} ==1 ));
%     [p,~]= signrank(gain_rt{r}.all(p_list_rt{r} ==1 ));
%     [p,~]= signrank(gain_r2{r}.all(stim_r2{r}.all ==1 ));
    title(['p = ',num2str(p,3)])
%     title([num2str(mean(gain_r2{r}.all(stim_r2{r}.all)),4),"   ",num2str(mean(gain_r2{r}.all(p_list_r2{r})),4)])
    
    % suppression
    subplot(2,1,2)
    edges = -2.0:0.1:2.0;
    histogram(gain_r2{r}.supp(logical(stim_r2{r}.supp) & abs(gain_r2{r}.supp)<1.5),edges,'facecolor','k')
    hold on
    histogram(gain_r2{r}.supp(p_list_supp{r}),edges,'facecolor',"g")
    
    if sum(logical(stim_r2{r}.supp) & abs(gain_r2{r}.all)<1.5) > 0
        xline(mean(gain_r2{r}.supp(logical(stim_r2{r}.supp) & abs(gain_r2{r}.supp)<1.5)),'-r')
        xline(median(gain_r2{r}.supp(logical(stim_r2{r}.supp) & abs(gain_r2{r}.supp)<1.5)),'--r')
    
    hold off
%     histogram(gain_rt{r}.all(p_list_rt{r}),edges,'facecolor',"#D95319")
%     hold on
%     histogram(gain_rt{r}.all(logical(not(p_list_rt{r}).*stim_rt{r}.all)),edges,'facecolor','k')
%     [p,~]= signrank(gain_r2{r}.supp(logical(stim_r2{r}.supp)&abs(gain_r2{r}.all)<1.5));
    [p,~]= signrank(gain_r2{r}.supp(p_list_supp{r} == 1));
    title(['p = ',num2str(p,3)])
    ylim([0,10])
    end
%     title(num2str(mean(gain_rt{r}.all(p_list_rt{r}))))
end


% mean(gain_r2{r}.all(stim_r2{r}.all))
% 
r = 1;
[T,p]= ranksum(gain_r2{r}.all(p_list_r2{r} ==1),gain_rt{r}.all(p_list_rt{r} ==1) )
% [p,~]= signrank(gain_r2{r}.all(logical(stim_r2{r}.supp)&abs(gain_r2{r}.all)<1.5))
% [p,~,stats]= signrank(gain_r2{r}.supp(stim_r2{r}.supp ==1))




% for r= 1:2
%     figure
%     edges = -1.5:0.05:1.5;
%     histogram(gain_r2{r}.supp(stim_r2{r}.supp ==1),edges,'facecolor','k')
%     hold on
%     histogram(gain_r2{r}.supp(logical((stim_r2{r}.supp ==1).*(pgain{r}.p < 0.05))),edges,'facecolor','g')
%     title(num2str(mean(gain_r2{r}.supp(stim_r2{r}.supp ==1)),3))
%     ylim([0,15])
% end
%% calculate nb of units encoding stim

five = unique([list_r1.five, list_r2.five, list_rt.five]);
ten = unique([list_r1.ten, list_r2.ten, list_rt.ten]);





%%
% pgain = ana_pregain(rate,Pool);
% gain_r2{1}.all = -gain_r2{1}.all;





gain_all_e ={};
gain_all_i = {};
for r = 1:2
    figure
    
    edges = -1:0.05:1;
    if r == 1
        gain_all_e{1} = -gain_all_e{1};
    end
    
    histogram(gain_all_e{r}(gain_all_e{r} ~=0),edges,'facecolor','k')
    hold on 
    
%     title([num2str(median(gain_all_e{r}(gain_all_e{r} ~=0)),4),"   ",num2str(mean(gain_all_e{r}(gain_all_e{r} ~=0)),4)])

end



for r = 1:2
[p,~] = signrank(gain_all_e{r}(gain_all_e{r} ~=0))
end

% sum((gain_all2{1} ~=0) | (gain_all2{2} ~=0))
%%
r =2
sum(stim_r2{1}.all) %.*stim_r2{2}.all)
sum(stim_r2{2}.all)
sum(stim_r2{1}.all.*stim_r2{2}.all)

sum(stim_r2{1}.supp) %.*stim_r2{2}.all)
sum(stim_r2{2}.supp)
sum(stim_r2{1}.supp.*stim_r2{2}.supp)



[gain_r2, stim_r2,rate,~] = ana_gain2(Pool,2,2);



%%


r = 1
median(gain_r2{r}.supp(stim_r2{r}.supp ==1))

for r= 1:2
    figure
    edges = -1.5:0.05:1.5;
    histogram(gain_r2{r}.all(stim_r2{r}.all),edges,'facecolor','k')
    hold on
    histogram(gain_r2{r}.all(logical(stim_r2{r}.all.*(pgain{r}.p < 0.05))),edges,'facecolor','r')
%     title(num2str(mean(gain_r2{r}.supp(stim_r2{r}.supp ==1)),3))
    ylim([0,10])
end

%%
% violinplot([gain_all{r}(1,:),gain_all{r}(2,:)],[ones(1,length(gain_all{r})),ones(1,length(gain_all{r}))*2], ...
%     'Orientation', 'vertical','ViolinColor' ,cmap{r},'ShowMean',true );
% ylim([-1.4,1.4])



% cmap = colormap(parula(70));
% c = zeros(length(Pool),3);
% for n = 1:length(Pool)
%     c(n,:) = cmap(Pool(n).ttr-200,:);
% end

% for gr = 1:max(ttr_list)

test = [];
for r = 1:2
    figure
    boxplot([gain_all{r}(2,:),gain_all{r}(1,:)],[ones(1,length(gain_all{r})),ones(1,length(gain_all{r}))*2])
    hold on
    
    for n = 1:length(gain_all{r})
%         if ttr_list(n_list{r}(n)) == gr
            if gain_all{r}(1,n) > 0
                plot([2;1],gain_all{r}([1:2],n),'-b') %'Color',c(n_list{r}(n),:))
                
%                 test = [test;[gain_all{r}(2,n),Pool(n_list{r}(n)).ttr]];
            elseif gain_all{r}(1,n) < 0
                plot([2;1],gain_all{r}([1:2],n),'-r') %'Color',c(n_list{r}(n),:),'LineStyle','--')
            end
%         end
    end
    ylim([-2.0,2.0])
%     colormap(parula(70))
%     colorbar
end
% end

% figure
% x = test(:,2)-200;
% y1 = test(:,1);
% mdl =fitlm(x,y1);
% scatter(x,y1);
% hold on 
% plot(mdl)
% title(mdl.Rsquared.Adjusted)


r=1;
test1 = gain_all{r}(1,(gain_all{r}(1,:) < 0));
test2 = gain_all{r}(2,(gain_all{r}(1,:) < 0));


% test1 = gain_r2{r}.all(p_list_r2{r});
% test2 = gain_rt{r}.all(p_list_rt{r});

[p, S] = signrank(test1,test2)



%%
a_list = {};
for r = 1:2
    a_list{r} = zeros(1,length(gain_all{r}));
    for n = 1:length(gain_all{r})
        if gain_all{r}(1,n) > 0
            %             if gain_all{r}(1,n)-gain_all{r}(2,n) > 0
            a_list{r}(1,n) = gain_all{r}(1,n)-gain_all{r}(2,n);
            %
            %                 a_list{r}(1,n) = 1;
            %             else
            %                 a_list{r}(1,n) = -1;
            %             end
        elseif gain_all{r}(1,n) <= 0
%             if gain_all{r}(1,n)-gain_all{r}(2,n) > 0
                a_list{r}(1,n) = -gain_all{r}(1,n)+gain_all{r}(2,n);
                %                 a_list{r}(1,n) = -1;
                %             else
                %                 a_list{r}(1,n) = 1;
%             end
            %         else
            %                 a_list{r}(1,n) = -1;
        end
        
    end
end

%%
r = 2;
ttr_n_list = [Pool(n_list{r}).ttr];
ttr_list = unique(ttr_n_list);
y = zeros(1,length(ttr_list));

figure
y_for_fit = [];
for tr = 1:length(ttr_list)
    scatter(ttr_n_list(1,(ttr_n_list == ttr_list(tr))),a_list{r}(1,(ttr_n_list == ttr_list(tr))));
    y(1,tr) = mean(a_list{r}(1,(ttr_n_list == ttr_list(tr))));
    y_for_fit = [y_for_fit,a_list{r}(1,(ttr_n_list == ttr_list(tr)))];
    hold on
end
x = ttr_n_list - 200;
plot(ttr_list,y)
ylim([-1,1])
figure
mdl =fitlm(x,a_list{r});
% scatter(x,y1);
% hold on 
plot(mdl)
title(mdl.Rsquared.Adjusted)

%% plot PSTH of selected units
s_ind =2;
if s_ind == 2
    phase = [2,4];
else
    phase = [1,3];
end
win = 10;
Y = {};
Y{1} = [];
Y{2} = [];
Y{3} = [];
n_list = find(not(stim_10{s_ind}.all >0)&stim_t{s_ind}.all >0);
% n_list = find((stim_20{s_ind}.all >0));

for n = n_list
    temp_max =  zeros(1,3);
    for p = 1:length(phase)
        temp_y = mean(rate.PSTH{n,phase(p)},1);
        temp_y = reshape(temp_y,win,[]);
        temp_y = mean(temp_y,1);
        Y{p} = [Y{p};temp_y];
        temp_max(p) = max(temp_y);
    end
    
    for p = 1:length(phase)
        Y{p}(end,:) = Y{p}(end,:)/max(temp_max);
    end
end





c = {'-g','-r','-b'};

sm = 10;

figure;
for p = 1:length(phase)
    y = mean(Y{p},1);
    ys = smoothdata(y,"gaussian",sm);
    err = std(Y{p},1) ...
        /sqrt(size(Y{p},1));
    err = smoothdata(err,"gaussian",sm);
    
    shadedErrorBar([1:length(ys)],ys,err,'lineProps',c{p})
    hold on
    xlim([100,400])
    ylim([0,0.4])
end


%% plot waveform of all n_list units
s_ind = 1;
figure
y_mean = [];
for n = 1:length(n_list{s_ind})
    if (gain_all{s_ind}(1,n) >0)
        y = mean(Pool(n_list{s_ind}(n)).waveforms,1)/max(mean(abs(Pool(n_list{s_ind}(n)).waveforms),1));
        plot(y-mean(y(1:10)))
        hold on
        y_mean = [y_mean;y-mean(y(1:10))];
    end
end
plot(nanmean(y_mean,1),'-k','Linewidth',2)



%% gain modulation, separate



[gain, stim] = ana_gain(Pool,2);

gain_real = gain;


listTV = {};

listTV{1} = load('D:\DATA\listTV2_AC_R1.mat');
listTV{2} = load('D:\DATA\listTV2_AC_RT.mat');


list_stim = {};
list_stim{1} = find(stim{1}.all == 1);
list_stim{2} = find(stim{2}.all == 1);

list_stim2 = {};
list_stim2{1} = unique([listTV{1}.five,listTV{2}.five])+1;
list_stim2{2} = unique([listTV{1}.ten,listTV{2}.ten])+1;

inter = {};
inter{1} = intersect(list_stim{1},list_stim2{1});
inter{2} = intersect(list_stim{2},list_stim2{2});


for n = 1:length(Pool)
    for r = 1:2
        if ~ismember(n,inter{r})
            gain{r}.pv(1,n) = 1;
        end
    end
end


n_list = [];
for n = list_stim{2}
    if ~ismember(n,inter{2})
        n_list = [n_list,n];
    end
end

% gain_rt = gain;
%% scatter plot RT vs R2   


for r = 1:2
    p_list_r2{r} = (0 < gain_r2{r}.pv) & (gain_r2{r}.pv < 0.05);
    p_list_rt{r} = (0 < gain_rt{r}.pv) & (gain_rt{r}.pv < 0.05);

end

cmap = colormap(jet(70));
c = zeros(length(Pool),1);
for n = 1:length(Pool)
    c(n) = Pool(n).ttr-200;
end

for r =1:2
    figure
    scatter(gain_r2{r}.all((p_list_r2{r}|p_list_rt{r})),gain_rt{r}.all((p_list_r2{r}|p_list_rt{r})), ...
        30,c((p_list_r2{r}|p_list_rt{r})),'filled')
    colormap(jet(70));
    ylim([-1.2,1.2])
    xlim([-1.2,1.2])
    colorbar
    hold on
    xline(0)
    yline(0)
    
end


%%

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
% 
% vd2 = [gain{1}.all( g_list{1}),gain{2}.all( g_list{2})];
% vc2 = [ones(sum(g_list{1}),1)*2;ones(sum(g_list{2}),1)*1];
% vc2 = vc2.';


color = {[0,0,1],[1,0,1],[0,0,1],[1,0,1]};
figure
% vs2 = violinplot(vd2,vc2,'Orientation', 'horizontal','ViolinColor',[0,0,0]);
hold on
vs1 = violinplot(vd,vc,'Orientation', 'horizontal','ViolinColor' ,[[1,0,1];[0,0,1]],'ShowMean',true );
hold off
% vs1(1,1).ScatterPlot.MarkerFaceColor = [1,1,1];

xlim([-1.2,1.2])