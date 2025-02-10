
function plot_raster(Pool, raster, rate, lick, lick_rate,fig_list)


stim_order = [1,2,7,8,3,4]; % R1_Go, R1_NoGo, RT_Go, RT_NoGo, R2_Go, R2_NoGo
c_list = {'r','b','r','b','r','b','r','b'};
ls_list = {'-','-','--','--',':',':',':',':'};
label = {'5_{r1}','10_{r1}','5_{r2}','10_{r2}','5_{CD}','10_{CD}','5_{CD}','10_{CD}'};

if isempty(fig_list)
    P_list = 1:length(Pool);
else
    P_list = fig_list;
end
    

for pp = 1:length(P_list)
    p = P_list(pp);
    
    
    f = figure(p);
    f.Position = [500 100 800 800];
    subplot(3,2,[1,3])
    stim_dur = 0.5;
    delay = 1.0;
    reward = 2.0;
    TotalReps= size(Pool(p).xb.trial_type,1);
    
%     nreps1 = cumsum([raster{p}.nreps(stim_order)]);
    for st = 1:length(raster{p})
        if isempty(raster{p}(st).nreps)
            raster{p}(st).nreps = 0;
        end
    end
    nreps1 = [raster{p}.nreps];
    nreps1 = cumsum(nreps1(stim_order));
    rectangle('Position',[0 0,stim_dur nreps1(end)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
    
    rectangle('Position',[(stim_dur + delay) 0,(stim_dur + delay+reward) nreps1(end)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
    hold on
%     nreps = nreps1([raster.stim{p}]);
%     yline(nreps1, color = 'r', linewidth = 2)
    r = 0;
    for st = stim_order
        plot(raster{p}(st).spikes,r+raster{p}(st).rep,'k.','MarkerSize',10);
        hold on
        r = r+ raster{p}(st).nreps;
    end
    hold off
    
    subplot(3,2,5)
    for st = stim_order
        plot(smoothdata(mean(rate.PSTH{p,st}),"gaussian",200), ...
                color = c_list{st}, linewidth = 2, linestyle = ls_list{st}, ...
                DisplayName = label{st})
        hold on
    end
    hold off
    
    subplot(3,2,[2,4])
        rectangle('Position',[0 0,stim_dur nreps1(end)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
    
    rectangle('Position',[(stim_dur + delay) 0,(stim_dur + delay+reward) nreps1(end)],'FaceColor',[0.9,0.9,0.9],'EdgeColor','none')
    hold on
%     nreps = nreps1([raster.stim{p}]);
    yline(nreps1, color = 'r', linewidth = 2)
    r = 0;
    for st = stim_order
        plot(lick{p}(st).licks,r+lick{p}(st).rep,'k.','MarkerSize',10);
        hold on
        r = r+ lick{p}(st).nreps;
    end
    
    
    hold off
    
        subplot(3,2,6)
    for st = stim_order
        plot(smoothdata(mean(lick_rate.PSTH{p,st}),"gaussian",200), ...
                color = c_list{st}, linewidth = 2, linestyle = ls_list{st}, ...
                DisplayName = label{st})
        hold on
    end
    hold off
    
    legend
    pause(0.1)
    sgtitle(['unit' num2str(p)])
    figname = fullfile('D:\DATA\figures',filesep, ['unit' num2str(p) '.png']);
    saveas(f,figname)
    close(f)
end



%% gather rate
% 
% R = {};
% for st = 1:4
%     R{st} = zeros(length(Pool),8000);
% 
%     for p  =1:length(Pool)
%         R{st}(p,:) = mean(rate.PSTH{p,st},1);
%     end
% end
% 
% 
% 
% figure
% pre= 2000;
% for st = 1:4
%     A = mean(R{st},1)-mean(R{st}(:,1:pre));
%     plot(smoothdata(A,"gaussian",200))
%     hold on
% end