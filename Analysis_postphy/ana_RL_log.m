

log  ={};

for i = 1:13
    for ii = 5:length(fromtxt{1})
    log{ii-4,i} = fromtxt{i}{ii};
        if i == 4
            log{ii-4,i} = floor(str2double(log{ii-4,i})/10);
        end
    end
end

%% 

trial_ind = find(strcmp([log{:,3}],"Sound"));
p_ind = find(strcmp(log(:,3),'Picture'));



t_ind = find([log{:,4}] == 31 | [log{:,4}] == 35 | [log{:,4}] == 41 | [log{:,4}] == 45);

h_ind = [log{t_ind,4}].';

hitcode= zeros(length(h_ind),4);
for tr  = 1:length(h_ind)
    if h_ind(tr) == 35
        hitcode(tr,1) = 1;
    elseif h_ind(tr) == 45
        hitcode(tr,2) = 1;
    elseif h_ind(tr) == 31
        hitcode(tr,3) = 1;
    elseif h_ind(tr) == 41
        hitcode(tr,4) = 1;
    end
end
    

t_ind = find([log{p_ind,4}]  < 100);
t_ind2 = find([log{p_ind,4}] > 100 & [log{p_ind,4}] < 300);
stim_idx = [log{trial_ind,4}];
t_idx = [log{p_ind(t_ind),4}];
t_idx2 = [log{p_ind(t_ind2),4}];

t_id = {};
t_id.left = t_idx(t_idx2 ==211);
t_id.right = t_idx(t_idx2 == 212);


Rate = {};
Rate.all = [];
Rate.right = [];
Rate.left = [];

h = 4; % Hit :1, Miss : 2, FA : 4
for t = 1:length(t_ind)-10
    Rate.all(t) = sum(t_idx(t:t+10)==h)/10;
end
for t = 1:length(t_id.left)-10
    Rate.left(t) = sum(t_id.left(t:t+10)==h)/10;
end
for t = 1:length(t_id.right)-10
    Rate.right(t) = sum(t_id.right(t:t+10)==h)/10;
end

plot(Rate.all)
hold on 
% plot(Rate.right)