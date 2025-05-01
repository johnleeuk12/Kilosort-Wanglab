%{
09/24/24 JHL
Analysis code and functions for gain-modulation index analysis in
manuscript. For figure 1d-h 
for details on how we calculate gain, read the function "ana_gain"
%} 

%% determining gain

Pool = load('Pool_IC.mat');
[gain_IC,stim_IC] = ana_gain(Pool);
Pool = load('Pool_AC.mat');   
[gain_AC,stim_AC] = ana_gain(Pool);
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


%% calculate mean gain for each area
sum((stim{1}.all+stim{2}.all)>0)
r = 1  
[~,p ] = kstest(gain_AC{r}.all(p_list2{r}))
mean(gain_AC{r}.all(p_list2{r}))



