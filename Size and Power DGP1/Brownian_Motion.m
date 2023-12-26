
n=5000;
R = 5000;
rv = nan(R,1);

pi0 = 0.25;

for r=1:R

e=randn(n,1);
mu = mean(e);
e_demean = e-mu;

den=(1/n)*sum((cumsum(e_demean/sqrt(n))).^2);
num = sum(e)/sqrt(n);

rv(r) = (num)/sqrt(den);


end

ptiles_rv = (prctile(rv,[1 2.5 5 10 90 95 97.5 99]))


