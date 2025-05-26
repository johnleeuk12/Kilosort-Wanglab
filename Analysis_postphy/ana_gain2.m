function [gain,stim,rate,raster2] = ana_gain2(Pool,rule,sd)
addpath('D:\GitHub\Kilosort-Wanglab\Analysis_postphy_core')
clear raster2 rate Lick2 Lick_rate
[raster2, rate, Lick2, Lick_rate] = gather_raster_ephys(2, 6, Pool);

plot_raster(Pool, raster2, rate, Lick2, Lick_rate,[20:30])
if rule == 2
    r_ind5 = 3;
    r_ind10 = 4;
elseif rule ==3
    r_ind5 = 7;
    r_ind10 = 8;
end

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
    stim{r}.supp = zeros(1,length(list_n));
    gain{r}.pv =  zeros(1,length(list_n));
    gain{r}.supp = zeros(1,length(list_n));
    supp{r}.sus =  zeros(1,length(list_n));
    supp{r}.all = zeros(1,length(list_n));
end
alpha = 2;
for n = 1:length(list_n)
    nn = list_n(n);
    for r = 1:2 % 5 or 10 kHz, analysed separately
        if r == 1
            y_go =  smoothdata(mean(rate.PSTH{nn,1},1),"gaussian",50);
            raw_go = rate.PSTH{nn,1};
            y_ng = smoothdata(mean(rate.PSTH{nn,r_ind5},1),"gaussian",50);
            raw_ng = rate.PSTH{nn,r_ind5};
        else
            y_go =  smoothdata(mean(rate.PSTH{nn,r_ind10},1),"gaussian",50);
            raw_go = rate.PSTH{nn,r_ind10};
            y_ng = smoothdata(mean(rate.PSTH{nn,2},1),"gaussian",50);
            raw_ng = rate.PSTH{nn,2};
            
        end
        d1 = y_go(2000:2100);
        d2 = y_go(500:1500);
        d3 = y_ng(2000:2100);
        d4 = y_ng(500:1500);
        
        bootstat = bootstrp(100,@mean,d1);
        bootstat2 = bootstrp(100,@mean,d2);
        [~,p1] = kstest2(bootstat,bootstat2);
        bootstat = bootstrp(100,@mean,d3);
        bootstat2 = bootstrp(100,@mean,d4);
        [~,p2] = kstest2(bootstat,bootstat2);        
        
        if p1 < 0.05 || p2 <0.05
            if mean(d1) > mean(d2) || mean(d3)>mean(d4)
                stim{r}.onset(n) = 1;
                gain{r}.onset(n) = ((mean(y_go(2000:2100))- mean(y_go(500:1500)))- ...
                    (mean(y_ng(2000:2100))- mean(y_ng(500:1500))))/ ...
                    (max([(mean(y_go(2000:2100))- mean(y_go(500:1500))), ...
                    (mean(y_ng(2000:2100))- mean(y_ng(500:1500)))]) + ...
                    alpha);
            end
        end
        
        d1 = y_go(2100:2450);
        d2 = y_go(500:1500);
        d3 = y_ng(2100:2450);
        d4 = y_ng(500:1500);
        
        bootstat = bootstrp(100,@mean,d1);
        bootstat2 = bootstrp(100,@mean,d2);
        [~,p1] = kstest2(bootstat,bootstat2);
        bootstat = bootstrp(100,@mean,d3);
        bootstat2 = bootstrp(100,@mean,d4);
        [~,p2] = kstest2(bootstat,bootstat2);  
        
        if p1 < 0.001 || p2 <0.001
            if mean(d1) > mean(d2) || mean(d3)>mean(d4)
                stim{r}.sus(n) = 1;
                gain{r}.sus(n) = ((mean(y_go(2000:2450))- mean(y_go(500:1500)))- ...
                    (mean(y_ng(2100:2450))- mean(y_ng(500:1500))))/ ...
                    (max([(mean(y_go(2100:2450))- mean(y_go(500:1500))), ...
                    (mean(y_ng(2100:2450))- mean(y_ng(500:1500)))]) + ...
                    alpha);
            elseif mean(d1) < mean(d2) && mean(d3) < mean(d4)
                if stim{r}.onset == 0
                    
                    stim{r}.supp(n) = 1;
                    
                    gain{r}.sus(n) = -((mean(y_go(2000:2450))- mean(y_go(500:1500)))- ...
                        (mean(y_ng(2100:2450))- mean(y_ng(500:1500))))/ ...
                        (max([(mean(y_go(2100:2450))- mean(y_go(500:1500))), ...
                        (mean(y_ng(2100:2450))- mean(y_ng(500:1500)))]) + ...
                        alpha);
                    gain{r}.supp(n) = gain{r}.sus(n);
                end
            end
        end

%        
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
%         if max([(mean(y_go(2100:2450))- mean(y_go(500:1500))), ...
%                 (mean(y_ng(2100:2450))- mean(y_ng(500:1500)))])...
%                 < 0
% %             supp{r}.all(n) = gain{r}.sus(n);
%             gain{r}.sus(n) = 0;
%         end
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
        
%         if median(y_go(2100:2450)) == 0 && median(y_ng(2100:2450)) == 0
%             gain{r}.sus(n) = 0;
%             stim{r}.sus(n) = 0;
%         end
        % removing units end
        
        if abs(gain{r}.onset(n)) > abs(gain{r}.sus(n))
            gain{r}.all(n) = gain{r}.onset(n);
            d = (mean(raw_go(:,2000:2100),2)-mean(y_go(500:1500)));
            d2 = (mean(raw_ng(:,2000:2100),2)-mean(y_ng(500:1500)));
            d  = d + 2*max(mean(y_go(500:1500)),mean(y_ng(500:1500)));
            d2 = d2 + 2*max(mean(y_go(500:1500)),mean(y_ng(500:1500)));
            
            p = spk_ttest2(d,d2);
            
            gain{r}.pv(n) = p;
        elseif abs(gain{r}.onset(n)) < abs(gain{r}.sus(n))
            gain{r}.all(n) = gain{r}.sus(n);
            d = (mean(raw_go(:,2100:2450),2)-mean(y_go(500:1500)));
            d2 = (mean(raw_ng(:,2100:2450),2)-mean(y_ng(500:1500)));
            d  = d + 2*max(mean(y_go(500:1500)),mean(y_ng(500:1500)));
            d2 = d2 + 2*max(mean(y_go(500:1500)),mean(y_ng(500:1500)));
            
            p = spk_ttest2(d,d2);
            gain{r}.pv(n) = p;
        end
        if stim{r}.supp(n) ==1
            d = (mean(raw_go(:,2100:2450),2)-mean(y_go(500:1500)));
            d2 = (mean(raw_ng(:,2100:2450),2)-mean(y_ng(500:1500)));
            d  = d + 2*max(mean(y_go(500:1500)),mean(y_ng(500:1500)));
            d2 = d2 + 2*max(mean(y_go(500:1500)),mean(y_ng(500:1500)));
            p = spk_ttest2(d,d2);
            gain{r}.pv_supp(n) = p;
        end
    end
end

for r = 1:2
    stim{r}.all = ((stim{r}.onset + stim{r}.sus)>0);
end