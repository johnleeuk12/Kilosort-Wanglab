function pgain = ana_pregain(rate,Pool)

N = length(Pool);
pgain = {};
for r =1:2
    pgain{r}.v = zeros(1,N);
    pgain{r}.p = zeros(1,N);
    
end

for nn = 1:N
    for r = 1:2
        if r ==1
            r1 = 1;
            r2 = 3;
        else
            r1 = 2;
            r2 = 4;
        end
        pre_r1 = mean(mean(rate.PSTH{nn,r1}(:,500:1500)));
        pre_r2 = mean(mean(rate.PSTH{nn,r2}(:,500:1500)));
        raw_r1 = smoothdata(rate.PSTH{nn,r1}(:,500:1500),2,"gaussian",5);
        raw_r2 = smoothdata(rate.PSTH{nn,r2}(:,500:1500),2,"gaussian",5);
        [h,p] = kstest2(reshape(raw_r1,[],1),reshape(raw_r2,[],1));
%         [p,~] = signrank(reshape(raw_r1,[],1),reshape(raw_r2,[],1));

        pgain{r}.p(nn) = p;
        pgain{r}.v(nn) = (pre_r2-pre_r1)/(max(pre_r2,pre_r1)+1);
    end
end


% for r = 1:2
%     figure
%     edges = -1:0.05:1;
%     histogram(pgain{r}.v,edges,'facecolor',"k")
%     hold on
%     %     histogram(gain_r2{r}.all(logical(not(p_list_r2{r}).*stim_r2{r}.all)),edges,'facecolor','k')
%     histogram(pgain{r}.v(pgain{r}.p < 0.01),edges,'facecolor','g')
% end
% 
% mean(pgain{1}.v)

% r = 2
% sum((pgain{r}.p < 0.05))
% sum((stim_r2{r}.supp ==1))
% sum(logical((stim_r2{r}.supp ==1).*(pgain{r}.p < 0.01)))
% 
% r = 2
% [p,t,stats] = signrank(pgain{r}.v)