


CR = zeros(2,length(Pool));
n_list = 1:length(Pool);
n_list = n_list(list3);
for n = n_list
    x = [C2_1.hit(:,n); C2_1.FA(:,n);C2_2.hit(:,n); C2_2.FA(:,n)];
    l = [L_1.hit(:,n); L_1.FA(:,n);L_2.hit(:,n); L_2.FA(:,n)];
    x = smoothdata(x,'gaussian',20);
    l = smoothdata(l,'gaussian',20);
    [CR(1,n),CR(2,n)] = corr(x,l);
end

CR_AC = CR(:,n_list);



CR_AC(1,(CR_AC(2,:)>0.05)) = 0;
CR_IC(1,(CR_IC(2,:)>0.05)) = 0;

randIC = rand(1,length(CR_IC))/4 - 0.125;
randAC = rand(1,length(CR_AC))/4 - 0.125;

figure
scatter(1+(ones(1,length(CR_AC))+randAC),abs(CR_AC(1,:)),10,'filled')
hold on
scatter(ones(1,length(CR_IC))+randIC,abs(CR_IC(1,:)),10,'filled')

boxplot(abs([CR_AC(1,:),CR_IC(1,:)]),[ones(1,length(CR_AC)),zeros(1,length(CR_IC))])

sum(CR_AC(2,:)>0.05)

[T,p] = kstest2(CR_AC(1,:),CR_IC(1,:))