function C3 = normC2(C2,Pool,spont) 

C2.hit = smoothdata(C2.hit,"gaussian",10); 
C2.FA = smoothdata(C2.FA,"gaussian",10);

for n = 1:length(Pool)
    C3.hit(:,n) = (C2.hit(:,n)-spont.mean(n,1))/(max(abs(C2.hit(:,n)-spont.mean(n,1)))+0.5);
    C3.FA(:,n) = (C2.FA(:,n)-spont.mean(n,1))/(max(abs(C2.FA(:,n)-spont.mean(n,1)))+0.5);
end
