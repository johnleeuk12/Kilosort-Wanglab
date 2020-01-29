%% Cf analysis 11/11/2019. 
% modify to include track information, and maybe plot onto a diagram?





%%
test1 = [M12E_unit_list.data{:,7}] ;
test1 = M12E_unit_list.data{:,7} ;
for i = 1:length(M12E_unit_list.data)
if strcmp(M12E_unit_list.data{:,7},'User: Trill_list.txt')
test = [test i]
end
end
Error using strcmp
Too many input arguments.
 test = [];
for i = 1:length(M12E_unit_list.data)
if strcmp(M12E_unit_list.data{i,7},'User: Trill_list.txt')
test = [test i];
end
end

for t = 1:length(test)
    M12E_unit_list.data{test(t),6} = 286;
end



%%
FR_norm = [];
% figure
%normalizing FR by spont rate
NN = size(SUrate_1,2);
for n = 1:NN
    FR_norm(n,:) = SUrate_1{n}.mean- SUrate_1{n}.spont; %./max(SUrate_1{n}.mean);
    
%     plot(FR_norm(n,:))
%     hold on
%     drawnow
end


% [Y,Sig,X] = svd(FR_norm,'econ');
% %                 sig = diag(Sig);%figure; semilogy(sig(sig>1),'kx-')
% k = 1:3;
% P = Y(:,k)*Sig(k,k)*X(:,k)';
% 
% for n = 1:53
%     plot(P(n,:))
%     hold on
%     pause
% end
% 
% data(ch,:) = mean(P,2).';

%find evoked responses for a certain range

% define stimuli range

% at least n_s amount of evoked responses


for i = 1:4
    xxx{i} = find(FR_norm(:,i+6)>0);
end

list = intersect(list, xxx{4});%, xxx{3}, xxx{4});
list = [xxx{1}; xxx{2}; xxx{3}; xxx{4}];
list = unique(list);



test = mean(FR_norm(list,:),1);
err = std(FR_norm,[],1);
err2 = std(FR_norm,[],2);
 [~,score,~] = pca(FR_norm);
scatter(score(:,1),score(:,2))
% x = 80:10:240;
SEM = err./sqrt(length(SUrate_1));
    ts = tinv([0.025  0.975],sqrt(length(SUrate_1)));  % T-Score
%     %     CI = mean2(output.rates_stim{f}) + ts(2)*SEM;                      % Confidence Intervals
    error = ts(2)*SEM;

errorbar(x,test,error);

%% 
% stim = 80:10:240;
% p =1;
N = length(rate_1.pre);

stim = 4:0.4:10;

FR_norm = [];
for n = 1:N
    counter = 0;
    std_spont = std(reshape(rate_1.pre{n},[],1),[],1)/sqrt(length((reshape(rate_1.pre{n},[],1))));
    ts = tinv([0.025  0.975],sqrt(length((reshape(rate_1.pre{n},[],1)))));
    SEM_spont = std_spont*ts(2);
    mean_spont = mean(reshape(rate_1.pre{n},[],1));
    
    %have at least one stim evoked response for a given stimuli
    for st = 1:length(stim)
        stim_rate = mean(rate_1.stim{n}(st,:),2);
        if abs(stim_rate) > mean_spont + SEM_spont*2
            counter =1;
        end
    end
    if counter == 1
        temp = mean(rate_1.stim{n},2).'-mean_spont;
        
        FR_norm = [FR_norm;temp ];
%         p = p+1;
    end
end

N = size(FR_norm,1);


FR_norm_2 = [];
stim_range = 3:12;
figure
for n = 1:N
    counter = 0;
    
    for st = 1:10
        if FR_norm(n,st+2)> 0
            counter = counter+1;
        end
    end
    if counter >2
        FR_norm_2 = [FR_norm_2; FR_norm(n,:)];
        plot(FR_norm(n,:))
        hold on
    end
end


SEM = std(FR_norm_2,[],1)/sqrt(length(FR_norm_2));
ts = tinv([0.025  0.975],sqrt(length(FR_norm_2)));

error = ts(2)*SEM;
errorbar(stim,mean(FR_norm_2),error)










% N_list_1 = unique([Pool_1.neuron_nb]);
% 
% for n = 1:length(N_list_1)
%     rec_list_1 = find([Pool_1.neuron_nb] ==N_list_1(n));
%     p = rec_list_1(1); %try end as well
%     nreps = size(rate_1.stim{1},2);
%     SUrate_1{n}.mean = mean(rate_1.stim{p},2); % -mean(rate_1.pre{p},2); %(Spikes/second)
%     SUrate_1{n}.error = std(rate_1.stim{p},1,2)/sqrt(nreps);
%     SUrate_1{n}.spont = mean(mean(rate_1.pre{p}));
%     StimDur = Pool_1(p).xb.stimulus_ch1(:,5)*1e-3;