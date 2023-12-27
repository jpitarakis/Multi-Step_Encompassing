function [T1,T1_alrv,T1_d,T1_d_alrv]=pred_encompass_dnorm(e1hat,e2hat,mu0)

% e1hat: n by 1 vector of out of sample forecast errors from model 1
% e2hat: n by 1 vector of out of sample forecast errors from model 2
% mu0: 


[n,~] = size(e1hat);
time_vec = (1:n)';
m0 = round(n*mu0);
w = (time_vec<=m0);

e1sq = e1hat.^2;
e2sq = e2hat.^2;
e12 = e1hat.*e2hat;
sigsq2 = mean(e2sq);
e2sq_demean = e2sq-sigsq2;


fm = ((1-2*mu0)^2)/(4*mu0*(1-mu0));
nw_lags = min(floor(1.2*n^(1/3)),n);

%d = e1sq - 0.5*((1/mu0)*e12.*w+(1/(1-mu0))*e12.*(1-w));
d = (e1sq-sigsq2) - 0.5*((1/mu0)*(e12-sigsq2).*w+(1/(1-mu0))*(e12-sigsq2).*(1-w));

sigsq_d = (d-mean(d))'*(d-mean(d))/n;
sigsq_d_nw = covnw(d,nw_lags,1);
sigsq_d_alrv = andrews_lrv(d);

phihatsq = (e2sq_demean)'*(e2sq_demean)/n;
phihatsq_nw = covnw(e2sq_demean,nw_lags,1);
phihatsq_alrv = andrews_lrv(e2sq_demean);

T1 = sqrt(n)*mean(d)/sqrt(fm*phihatsq);

T1_alrv = sqrt(n)*mean(d)/sqrt(fm*phihatsq_alrv);

T1_d = sqrt(n)*mean(d)/sqrt(sigsq_d);

T1_d_alrv = sqrt(n)*mean(d)/sqrt(sigsq_d_alrv);


