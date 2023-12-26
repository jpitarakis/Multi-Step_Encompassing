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

T1 = nan(R,sm,sh,sk);
T1_nw=nan(R,sm,sh,sk);
T1_alrv=nan(R,sm,sh,sk);

T1_d = nan(R,sm,sh,sk);
T1_d_nw=nan(R,sm,sh,sk);
T1_d_alrv=nan(R,sm,sh,sk);

size_T1_10 = nan(1,sm,sh,sk);
size_T1_nw_10 = nan(1,sm,sh,sk);
size_T1_alrv_10 = nan(1,sm,sh,sk);

size_T1_d_10 = nan(1,sm,sh,sk);
size_T1_d_nw_10 = nan(1,sm,sh,sk);
size_T1_d_alrv_10 = nan(1,sm,sh,sk);

theta = 0.5;
rng = 171966;

for k = 1:2
    N = N_vec(k);
    T = T_vec(k);
    
%rand('twister',312)
alpha = 0.2+rand(r,1)*0.6;
rho = 0.3+rand(N,1)*0.5;
%seed=123456;
%stream = RandStream('mrg32k3a','Seed',seed);


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

Y(h+1:n)=1.25+0.5*Y(1:n-h)+0*f_reg(1:n-h)+v(h+1:n);


% Estimating the factor (assume r=1 is known)

Xs = standard(X);
XX=Xs*Xs';
[Fhat0,eigval,Fhat1]=svd(XX');
Fhat = Fhat0(:,1:1)*sqrt(T+h+lags);
fhat_reg = Fhat0(:,1);


ehat1 = recursive_hstep_fast(Y(h+1:n),[Y(h+1:n)],0.25,h);
ehat2 = recursive_hstep_fast(Y(h+1:n),[Y(h+1:n),fhat_reg(h+1:n)],0.25,h);


for m=1:sm
[T1(b,m,d,k),T1_nw(b,m,d,k),T1_alrv(b,m,d,k),...
    T1_d(b,m,d,k),T1_d_nw(b,m,d,k),T1_d_alrv(b,m,d,k)]=pred_encompass_dnorm(ehat1,ehat2,mu_vec(m));
end

end

end
end

for k = 1:sk
for d = 1:sh
for m = 1:sm
            size_T1_10(1,m,d,k) = sum(T1(:,m,d,k)>1.2816)/R;
            size_T1_nw_10(1,m,d,k) = sum(T1_nw(:,m,d,k)>1.2816)/R;
            size_T1_alrv_10(1,m,d,k) = sum(T1_alrv(:,m,d,k)>1.2816)/R;
            size_T1_d_10(1,m,d,k) = sum(T1_d(:,m,d,k)>1.2816)/R;
            size_T1_d_nw_10(1,m,d,k) = sum(T1_d_nw(:,m,d,k)>1.2816)/R;
            size_T1_d_alrv_10(1,m,d,k) = sum(T1_d_alrv(:,m,d,k)>1.2816)/R;
end
end
end


% d based outcomes

mat1a_d_T10 = ["$\mu_{0}$","0.30","0.35","0.40","0.45"];
mat2a_d_T10 = ["","","(N,T)=(100,250)","",""];
mat3a_d_T10 = ["h=1",size_T1_d_10(:,:,1,1);"h=4",size_T1_d_10(:,:,2,1);"h=12", size_T1_d_10(:,:,3,1);"h=24", size_T1_d_10(:,:,4,1)];
mat4a_d_T10 = ["","","(N,T)=(500,500)","",""];
mat5a_d_T10 = ["h=1",size_T1_d_10(:,:,1,2);"h=4",size_T1_d_10(:,:,2,2);"h=12", size_T1_d_10(:,:,3,2);"h=24", size_T1_d_10(:,:,4,2)];

T10_d_all = [mat1a_d_T10;mat2a_d_T10;mat3a_d_T10;mat4a_d_T10;mat5a_d_T10];


% nw

mat1a_d_nw_T10 = ["$\mu_{0}$","0.30","0.35","0.40","0.45"];
mat2a_d_nw_T10 = ["","","(N,T)=(100,250)","",""];
mat3a_d_nw_T10 = ["h=1",size_T1_d_nw_10(:,:,1,1);"h=4",size_T1_d_nw_10(:,:,2,1);"h=12", size_T1_d_nw_10(:,:,3,1);"h=24", size_T1_d_nw_10(:,:,4,1)];
mat4a_d_nw_T10 = ["","","(N,T)=(500,500)","",""];
mat5a_d_nw_T10 = ["h=1",size_T1_d_nw_10(:,:,1,2);"h=4",size_T1_d_nw_10(:,:,2,2);"h=12", size_T1_d_nw_10(:,:,3,2);"h=24", size_T1_d_nw_10(:,:,4,2)];

T10_d_nw_all = [mat1a_d_nw_T10;mat2a_d_nw_T10;mat3a_d_nw_T10;mat4a_d_nw_T10;mat5a_d_nw_T10]

% alrv

mat1a_d_alrv_T10 = ["$\mu_{0}$","0.30","0.35","0.40","0.45"];
mat2a_d_alrv_T10 = ["","","(N,T)=(100,250)","",""];
mat3a_d_alrv_T10 = ["h=1",size_T1_d_alrv_10(:,:,1,1);"h=4",size_T1_d_alrv_10(:,:,2,1);"h=12", size_T1_d_alrv_10(:,:,3,1);"h=24", size_T1_d_alrv_10(:,:,4,1)];
mat4a_d_alrv_T10 = ["","","(N,T)=(500,500)","",""];
mat5a_d_alrv_T10 = ["h=1",size_T1_d_alrv_10(:,:,1,2);"h=4",size_T1_d_alrv_10(:,:,2,2);"h=12", size_T1_d_alrv_10(:,:,3,2);"h=24", size_T1_d_alrv_10(:,:,4,2)];

T10_d_alrv_all = [mat1a_d_alrv_T10;mat2a_d_alrv_T10;mat3a_d_alrv_T10;mat4a_d_alrv_T10;mat5a_d_alrv_T10];


%xlswrite('Montec_DGP2_size.xlsx',T10_d_all,'dgp2-d-T10','A1');
%xlswrite('Montec_DGP2_size.xlsx',T10_d_nw_all,'dgp2-d-T10nw','A1');
%xlswrite('Montec_DGP2_size.xlsx',T10_d_alrv_all,'dgp2-d-T10alrv','A1');


