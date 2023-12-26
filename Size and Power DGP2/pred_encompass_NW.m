function [T1,T1_nw]=pred_encompass_NW(e1hat,e2hat,mu0)

% e1hat: n by 1 vector of out of sample forecast errors from model 1
% e2hat: n by 1 vector of out of sample forecast errors from model 2
% mu0: 


[n,~] = size(e1hat);
time_vec = (1:n)';

e1sq = e1hat.^2;
e2sq = e2hat.^2;
e12 = e1hat.*e2hat;
e2sq_demean = e2sq-mean(e2sq);
m0 = round(n*mu0);

w = (time_vec<=m0);

fm = ((1-2*mu0)^2)/(4*mu0*(1-mu0));
nw_lags = min(floor(1.2*n^(1/3)),n);


d = e1sq - 0.5*((1/mu0)*e12.*w+(1/(1-mu0))*e12.*(1-w));

% T1, T1_nw (robust) 

phihatsq = (e2sq_demean)'*(e2sq_demean)/n;
phihatsq_nw = covnw(e2sq_demean,nw_lags,1);

T1 = sqrt(n)*mean(d)/sqrt(fm*phihatsq);
T1_nw = sqrt(n)*mean(d)/sqrt(fm*phihatsq_nw);






