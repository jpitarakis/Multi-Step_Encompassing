clear;
clc;

%tic
R=10000;

%N = 100;
%T=250;

N_vec=[100,500];
T_vec=[250;500];

[~,sk]=size(N_vec);

r = 1; % number of factors
lags = 1;
sigv=1;
mu_vec = [0.30;0.35;0.40;0.45];
[sm,~]=size(mu_vec);
h_vec = [1;4;12;24];
[sh,~]=size(h_vec);

beta2_vec = [0.00;0.10;0.20;0.30;0.35;0.40;0.45;0.50;0.55;0.60];
[sb,~] = size(beta2_vec);

T1 = nan(R,sm,sh,sb,sk);
T1_nw=nan(R,sm,sh,sb,sk);
T1_alrv=nan(R,sm,sh,sb,sk);

T1_d = nan(R,sm,sh,sb,sk);
T1_d_nw=nan(R,sm,sh,sb,sk);
T1_d_alrv=nan(R,sm,sh,sb,sk);

size_T1_10 = nan(sb,sm,sh,sk);

%size_T1_10 = nan(1,sm,sh,sb,sk);
size_T1_nw_10 = nan(1,sm,sh,sb,sk);
size_T1_alrv_10 = nan(1,sm,sh,sb,sk);

size_T1_d_10 = nan(1,sm,sh,sb,sk);
size_T1_d_nw_10 = nan(1,sm,sh,sb,sk);
size_T1_d_alrv_10 = nan(1,sm,sh,sb,sk);

theta = 0.5;
rng = 171966;

for k = 1:2
    N = N_vec(k);
    T = T_vec(k);
    
alpha = 0.2+rand(r,1)*0.6;
rho = 0.3+rand(N,1)*0.5;


for g = 1:sb
beta2 = beta2_vec(g);

for d = 1:sh
    
h = h_vec(d);
MAparams = theta.^((1:h-1)');
n = T+h+lags;

parfor b=1:R
 
eps = randn(T+h+lags,1);   
v = armaxfilter_simulate(eps,0,[],[],h-1,MAparams); 
u = randn(T+h+lags,r); 
F = u; % T x r

for j = 1:r
    F(:,j) = filter(1,[1,-alpha(j)],u(:,j)); %AR(1) process
end

eps = randn(T+h+lags,N);
e = eps;

for j=1:N
    e(:,j)=filter(1,[1,-rho(j)],eps(:,j)); 
end

lambda = randn(N,r)*sqrt(r); %factor lodadings
X = F*lambda'+e*sqrt(r); %generate panel data

f_reg = F(:,1);
[n1,~]=size(f_reg);

Y=zeros(n,1);

Y(h+1:n)=1.25+0.5*Y(1:n-h)+beta2*f_reg(1:n-h)+v(h+1:n);


% Estimating the factor (assume r=1 is known)

Xs = standard(X);
XX=Xs*Xs';
[Fhat0,eigval,Fhat1]=svd(XX');
Fhat = Fhat0(:,1:1)*sqrt(T+h+lags);
fhat_reg = Fhat0(:,1);


ehat1 = recursive_hstep_fast(Y(h+1:n),[Y(h+1:n)],0.25,h);
ehat2 = recursive_hstep_fast(Y(h+1:n),[Y(h+1:n),fhat_reg(h+1:n)],0.25,h);


for m=1:sm
[T1(b,m,d,g,k),T1_nw(b,m,d,g,k),T1_alrv(b,m,d,g,k),...
    T1_d(b,m,d,g,k),T1_d_nw(b,m,d,g,k),T1_d_alrv(b,m,d,g,k)]=pred_encompass_dnorm(ehat1,ehat2,mu_vec(m));
end

end

end
end

end



for g = 1:sb
for k = 1:sk
for d = 1:sh
for m = 1:sm
           size_T1_10(g,m,d,k) = sum(T1(:,m,d,g,k)>1.2816)/R;
           size_T1_nw_10(g,m,d,k) = sum(T1_nw(:,m,d,g,k)>1.2816)/R;
           size_T1_alrv_10(g,m,d,k) = sum(T1_alrv(:,m,d,g,k)>1.2816)/R;
          size_T1_d_10(g,m,d,k) = sum(T1_d(:,m,d,g,k)>1.2816)/R;
          size_T1_d_nw_10(g,m,d,k) = sum(T1_d_nw(:,m,d,g,k)>1.2816)/R;
          size_T1_d_alrv_10(g,m,d,k) = sum(T1_d_alrv(:,m,d,g,k)>1.2816)/R;
end
end
end
end

str_beta = ["$\beta_{2}$=0.0";"$\beta_{2}$=0.10";"$\beta_{2}$=0.20";"$\beta_{2}$=0.30";"$\beta_{2}$=0.35";"$\beta_{2}$=0.40";"$\beta_{2}$=0.45";"$\beta_{2}$=0.50";"$\beta_{2}$=0.55";"$\beta_{2}$=0.60";];
str_h1 = ["","","h=1","",""];
str_h4 = ["","","h=4","",""];
str_h12 = ["","","h=12","",""];
str_h24 = ["","","h=24","",""];
str_nt1 = ["","","$(N,T)=(100,250)$","",""];
str_nt2 = ["","","$(N,T)=(500,500)$","",""];
str_mu = ["$\mu_{0}","0.30","0.35","0.40","0.45"];

mat_T10_k1 = [str_nt1;str_mu;str_h1;str_beta,size_T1_10(:,:,1,1);str_h4;str_beta,size_T1_10(:,:,2,1);str_h12;str_beta,size_T1_10(:,:,3,1);str_h24;str_beta,size_T1_10(:,:,4,1)];
mat_T10_k2 = [str_nt2;str_mu;str_h1;str_beta,size_T1_10(:,:,1,2);str_h4;str_beta,size_T1_10(:,:,2,2);str_h12;str_beta,size_T1_10(:,:,3,2);str_h24;str_beta,size_T1_10(:,:,4,2)];
bigmat_T10 = [mat_T10_k1,mat_T10_k2];

mat_T10_nw_k1 = [str_nt1;str_mu;str_h1;str_beta,size_T1_nw_10(:,:,1,1);str_h4;str_beta,size_T1_nw_10(:,:,2,1);str_h12;str_beta,size_T1_nw_10(:,:,3,1);str_h24;str_beta,size_T1_nw_10(:,:,4,1)];
mat_T10_nw_k2 = [str_nt2;str_mu;str_h1;str_beta,size_T1_nw_10(:,:,1,2);str_h4;str_beta,size_T1_nw_10(:,:,2,2);str_h12;str_beta,size_T1_nw_10(:,:,3,2);str_h24;str_beta,size_T1_nw_10(:,:,4,2)];
bigmat_T10_nw = [mat_T10_nw_k1,mat_T10_nw_k2];

mat_T10_alrv_k1 = [str_nt1;str_mu;str_h1;str_beta,size_T1_alrv_10(:,:,1,1);str_h4;str_beta,size_T1_alrv_10(:,:,2,1);str_h12;str_beta,size_T1_alrv_10(:,:,3,1);str_h24;str_beta,size_T1_alrv_10(:,:,4,1)];
mat_T10_alrv_k2 = [str_nt2;str_mu;str_h1;str_beta,size_T1_alrv_10(:,:,1,2);str_h4;str_beta,size_T1_alrv_10(:,:,2,2);str_h12;str_beta,size_T1_alrv_10(:,:,3,2);str_h24;str_beta,size_T1_alrv_10(:,:,4,2)];
bigmat_T10_alrv = [mat_T10_alrv_k1,mat_T10_alrv_k2];

% d versions

mat_T10_d_k1 = [str_nt1;str_mu;str_h1;str_beta,size_T1_d_10(:,:,1,1);str_h4;str_beta,size_T1_d_10(:,:,2,1);str_h12;str_beta,size_T1_d_10(:,:,3,1);str_h24;str_beta,size_T1_d_10(:,:,4,1)];
mat_T10_d_k2 = [str_nt2;str_mu;str_h1;str_beta,size_T1_d_10(:,:,1,2);str_h4;str_beta,size_T1_d_10(:,:,2,2);str_h12;str_beta,size_T1_d_10(:,:,3,2);str_h24;str_beta,size_T1_d_10(:,:,4,2)];
bigmat_d_T10 = [mat_T10_d_k1,mat_T10_d_k2];

mat_T10_d_nw_k1 = [str_nt1;str_mu;str_h1;str_beta,size_T1_d_nw_10(:,:,1,1);str_h4;str_beta,size_T1_d_nw_10(:,:,2,1);str_h12;str_beta,size_T1_d_nw_10(:,:,3,1);str_h24;str_beta,size_T1_d_nw_10(:,:,4,1)];
mat_T10_d_nw_k2 = [str_nt2;str_mu;str_h1;str_beta,size_T1_d_nw_10(:,:,1,2);str_h4;str_beta,size_T1_d_nw_10(:,:,2,2);str_h12;str_beta,size_T1_d_nw_10(:,:,3,2);str_h24;str_beta,size_T1_d_nw_10(:,:,4,2)];
bigmat_T10_d_nw = [mat_T10_d_nw_k1,mat_T10_d_nw_k2];


mat_T10_d_alrv_k1 = [str_nt1;str_mu;str_h1;str_beta,size_T1_d_alrv_10(:,:,1,1);str_h4;str_beta,size_T1_d_alrv_10(:,:,2,1);str_h12;str_beta,size_T1_d_alrv_10(:,:,3,1);str_h24;str_beta,size_T1_d_alrv_10(:,:,4,1)];
mat_T10_d_alrv_k2 = [str_nt2;str_mu;str_h1;str_beta,size_T1_d_alrv_10(:,:,1,2);str_h4;str_beta,size_T1_d_alrv_10(:,:,2,2);str_h12;str_beta,size_T1_d_alrv_10(:,:,3,2);str_h24;str_beta,size_T1_d_alrv_10(:,:,4,2)];
bigmat_T10_d_alrv = [mat_T10_d_alrv_k1,mat_T10_d_alrv_k2];

%xlswrite('Montec_DGP2_power.xlsx',bigmat_T10,'dgp2-T10-power','A1');
%xlswrite('Montec_DGP2_power.xlsx',bigmat_T10_nw,'dgp2-T10-nw-power','A1');
%xlswrite('Montec_DGP2_power.xlsx',bigmat_T10_alrv,'dgp2-T10-alrv-power','A1');
%xlswrite('Montec_DGP2_power.xlsx',bigmat_d_T10,'dgp2-T10-d-power','A1');
%xlswrite('Montec_DGP2_power.xlsx',bigmat_T10_d_nw,'dgp2-T10-nw-dpower','A1');
%xlswrite('Montec_DGP2_power.xlsx',bigmat_T10_d_alrv,'dgp2-T10-alrv-d-power','A1');




