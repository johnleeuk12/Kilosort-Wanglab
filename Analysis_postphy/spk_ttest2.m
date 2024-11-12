function p = spk_ttest2(d,d2)


pd = fitdist(d,'normal');
pd2 = fitdist(d2,'normal');
t = abs((pd.mu-pd2.mu))/sqrt(pd.sigma^2/length(d)+pd2.sigma^2/length(d2));
df = (pd.sigma^2/length(d)+pd2.sigma^2/length(d2))^2/...
    ((pd.sigma^2/length(d))^2/(length(d)-1)+(pd2.sigma^2/length(d2))^2/(length(d2)-1));
p = 1-tcdf(t,df);