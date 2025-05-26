function [gain,stim,rate,raster2] = ana_gain3(Pool,rule,sd,sd2,list_n,list_np,list_nn)
addpath('D:\GitHub\Kilosort-Wanglab\Analysis_postphy_core')
clear raster2 rate Lick2 Lick_rate
[raster2, rate, Lick2, Lick_rate] = gather_raster_ephys(2, 6, Pool,list_n);

% plot_raster(Pool, raster2, rate, Lick2, Lick_rate,[20:30])
if rule == 2
    r_ind5 = 3;
    r_ind10 = 4;
elseif rule ==3
    r_ind5 = 7;
    r_ind10 = 8;
end

if isempty(list_n)
    list_n = 1:length(Pool);
end
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
    gain{r}.pv_supp = zeros(1,length(list_n));
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
        
        if ismember(nn,list_np{1,r})||ismember(nn,list_np{rule,r})
%         if (mean(y_go(2000:2100)) > mean(y_go(500:1500))+sd*std(mean(raw_go(:,500:1500),2))) ...
%                 || (mean(y_ng(2000:2100)) > mean(y_ng(500:1500))+sd*std(mean(raw_ng(:,500:1500),2)))
            if mean(y_go(2000:2100)) > mean(y_go(500:1500)) || mean(y_ng(2000:2100)) > mean(y_ng(500:1500))
                stim{r}.onset(n) = 1;
    %             stim{r}.all(n) = 1;
                gain{r}.onset(n) = ((mean(y_go(2000:2100))- mean(y_go(500:1500)))- ...
                    (mean(y_ng(2000:2100))- mean(y_ng(500:1500))))/ ...
                    (max([(mean(y_go(2000:2100))- mean(y_go(500:1500))), ...
                    (mean(y_ng(2000:2100))- mean(y_ng(500:1500)))]) + ...
                    alpha);
            end
%         end
%         
%         if (mean(y_go(2100:2450)) > mean(y_go(500:1500))+sd*std(mean(raw_go(:,500:1500),2))) ...
%                 || (mean(y_ng(2100:2450)) > mean(y_ng(500:1500))+sd*std(mean(raw_ng(:,500:1500),2)))
            if mean(y_go(2100:2450)) > mean(y_go(500:1500))|| mean(y_ng(2100:2450)) > mean(y_ng(500:1500))
                stim{r}.sus(n) = 1;
                gain{r}.sus(n) = ((mean(y_go(2000:2450))- mean(y_go(500:1500)))- ...
                    (mean(y_ng(2100:2450))- mean(y_ng(500:1500))))/ ...
                    (max([(mean(y_go(2100:2450))- mean(y_go(500:1500))), ...
                    (mean(y_ng(2100:2450))- mean(y_ng(500:1500)))]) + ...
                    alpha);            
%             if (mean(y_go(2000:2450))- mean(y_go(500:1500))) < 0 && ((mean(y_ng(2000:2450))- mean(y_ng(500:1500))) <0)
%                 gain{r}.sus(n) = -gain{r}.sus(n);
%             end
            end
        end
%         
        if ismember(nn,list_nn{1})|| ismember(nn,list_nn{rule})
            if stim{r}.onset(n) == 0
%             if (mean(y_go(2100:2450)) < mean(y_go(500:1500))-sd2*std(mean(raw_go(:,500:1500),2))) ...
%                     && (mean(y_ng(2100:2450)) < mean(y_ng(500:1500))-sd2*std(mean(raw_ng(:,500:1500),2)))
                %             stim{r}.sus(n) = 1;
                stim{r}.supp(n) = 1;
                gain{r}.supp(n) = -((mean(y_go(2000:2450))- mean(y_go(500:1500)))- ...
                    (mean(y_ng(2100:2450))- mean(y_ng(500:1500))))/ ...
                    (min([(mean(y_go(2100:2450))- mean(y_go(500:1500))), ...
                    (mean(y_ng(2100:2450))- mean(y_ng(500:1500)))]) - ...
                    alpha);
%             end

%             if (mean(y_go(2100:2450)) < mean(y_go(500:1500))-sd2*std(mean(raw_go(:,500:1500),2))) ...
%                     && not(mean(y_ng(2100:2450)) > mean(y_ng(500:1500))+sd*std(mean(raw_ng(:,500:1500),2)))
%                 %             stim{r}.sus(n) = 1;
%                 stim{r}.supp(n) = 1;
%                 gain{r}.supp(n) = -((mean(y_go(2000:2450))- mean(y_go(500:1500)))- ...
%                     (mean(y_ng(2100:2450))- mean(y_ng(500:1500))))/ ...
%                     (min([(mean(y_go(2100:2450))- mean(y_go(500:1500))), ...
%                     (mean(y_ng(2100:2450))- mean(y_ng(500:1500)))]) - ...
%                     alpha);
%             end
% 
% 
%             if not(mean(y_go(2100:2450)) > mean(y_go(500:1500))+sd*std(mean(raw_go(:,500:1500),2))) ...
%                     && (mean(y_ng(2100:2450)) < mean(y_ng(500:1500))-sd2*std(mean(raw_ng(:,500:1500),2)))
%                 %             stim{r}.sus(n) = 1;
%                 stim{r}.supp(n) = 1;
%                 gain{r}.supp(n) = -((mean(y_go(2000:2450))- mean(y_go(500:1500)))- ...
%                     (mean(y_ng(2100:2450))- mean(y_ng(500:1500))))/ ...
%                     (min([(mean(y_go(2100:2450))- mean(y_go(500:1500))), ...
%                     (mean(y_ng(2100:2450))- mean(y_ng(500:1500)))]) - ...
%                     alpha);
%             end
            end
        end

        
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