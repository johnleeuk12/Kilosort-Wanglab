


st_list = zeros(4, length(Pool));
Pre = 2*1e3;
stimdur = 0.5*1e3;

for n  = 1:length(Pool)
    for st = 1:4
        spont = rate.PSTH{st,n}(:,0:Pre);
        stim = rate.PSTH{st,n}(Pre:Pre+stimdur); 
    
    
    end
end
    
C4.IC = get_normC2(C2,spont,list_pre,'IC','FA');
C4.AC = get_normC2(C2,spont,list_pre,'AC','FA');





bin = 5;
edges = -1e3:bin:1e3;

figure
ys = [];
x = edges(1:end-1).';
for i = 1:length(mean(C4.IC,1))
    ys(i) = gaussian_kern_reg(x(i),x,mean(C4.IC,1).',20);
end
sd= std(C4.IC,[],1)/sqrt(size(C4.IC,1));
shadedErrorBar(edges(1:end-1),ys,sd,'lineProps','r')
% plot(mean(C4.AC,1))
hold on
% plot(mean(C4.IC,1))
ys = [];
x = edges(1:end-1).';
for i = 1:length(mean(C4.AC,1))
    ys(i) = gaussian_kern_reg(x(i),x,mean(C4.AC,1).',20);
end
sd= std(C4.AC,[],1)/sqrt(size(C4.AC,1));
shadedErrorBar(edges(1:end-1),ys,sd,'lineProps','b')


% figure
% plot(mean(C2.AC.hit(:,list_pre.AC(:,1)),2)) 
% hold on
% plot(mean(C2.IC.hit(:,list_pre.IC(:,1)),2))
% yline(mean(spont.IC.mean(:,1)),'--r')
% yline(mean(spont.AC.mean(:,1)),'-b')

