function C4 = get_normC2(C2,spont, list_pre,a_ind,b_ind)
% bin = 5;
% edges = -1e3:bin:1e3;

switch a_ind
    case 'AC'
        C2 = C2.AC;
        list_pre = list_pre.AC;
        spont = spont.AC;
    case 'IC'
        C2 = C2.IC;
        list_pre = list_pre.IC;
        spont = spont.IC;
end

switch b_ind
    case 'Hit'
        list_pre = list_pre(:,1);
        C2 = C2.hit;
    case 'FA'
        list_pre = list_pre(:,2);
        C2 = C2.FA;
end


C4 = C2(:,list_pre).'-spont.mean(list_pre,1);
% C4 =movmean(C4,5,2,'Endpoints','shrink');

% C4 = C2(:,list_pre)-mean(C2(1:25,list_pre),1);
% C4 =movmean(C4,5,1,'Endpoints','shrink');
% C4 =C4.';


alpha =1;
for n = 1:size(C4,1)
    C4(n,:) = C4(n,:)/(max(C4(n,:))+alpha);
end