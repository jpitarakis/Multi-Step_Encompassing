
R=10000;
pi0=0.25;

% dgp 

n = 500;
k0 = round(n*pi0);

mu_w = zeros(2,1);
Sigma_w = diag([1;0.25]);

mu_vec = [0.30;0.35;0.40;0.45];
rho_vec = [0.25;0.90;0.95];
h_vec = [1;4;12;24];

beta1 = 0.30;
beta2_vec = [0.00;0.10;0.20;0.30;0.35;0.40;0.45;0.50;0.55;0.60];

[sn,~] = size(beta2_vec);
[sm,~]=size(mu_vec);
[sr,~]=size(rho_vec);
[sh,~]=size(h_vec);

theta = 0.5;

% Initialisations

T1 = nan(R,sm,sr,sh,sn);
T1_nw=nan(R,sm,sr,sh,sn);
T1_alrv=nan(R,sm,sr,sh,sn);

T1_d = nan(R,sm,sr,sh,sn);
T1_d_nw=nan(R,sm,sr,sh,sn);
T1_d_alrv=nan(R,sm,sr,sh,sn);

size_T1_10 = nan(sn,sm,sr,sh);
size_T1_nw_10 = nan(sn,sm,sr,sh);
size_T1_alrv_10 = nan(sn,sm,sr,sh);

size_T1_d_10 = nan(sn,sm,sr,sh);
size_T1_nw_d_10 = nan(sn,sm,sr,sh);
size_T1_alrv_d_10 = nan(sn,sm,sr,sh);

tic

rng = 171966;


for l = 1:sn

beta2 = beta2_vec(l);
x = zeros(n,1);
y = zeros(n,1);


for s = 1:sr
rho = rho_vec(s);


for d = 1:sh
h = h_vec(d);
MAparams = theta.^((1:h-1)');

for j=1:R

w = mvnrnd(mu_w,Sigma_w,n);
eps = w(:,1);
v = w(:,2);

u = armaxfilter_simulate(eps,0,[],[],h-1,MAparams); 

for i=h+1:n   
    x(i) = rho*x(i-1)+v(i);
    y(i) = beta1*y(i-h)+beta2*x(i-h)+u(i);
end

data = [y(h+1:n),y(h+1:n),x(h+1:n)];

y = data(:,1);
X = data(:,2:3);


ehat1=recursive_hstep_fast(y,X(:,1),pi0,h);
ehat2=recursive_hstep_fast(y,X,pi0,h);

for m=1:sm
[T1(j,m,s,d,l),T1_nw(j,m,s,d,l),T1_alrv(j,m,s,d,l),...
    T1_d(j,m,s,d,l),T1_d_nw(j,m,s,d,l),T1_d_alrv(j,m,s,d,l)]=pred_encompass_dnorm(ehat1,ehat2,mu_vec(m));
end

end

end
end
end

toc
quant = [0.01, 0.025, 0.05, 0.10, 0.90, 0.95, 0.975, 0.99];
norm_cvs = norminv(quant);

for m = 1:sm
    for l = 1:sn
        for s = 1:sr
            for d  = 1:sh
            size_T1_10(l,m,s,d) = sum(T1(:,m,s,d,l)>1.2816)/R;
            size_T1_nw_10(l,m,s,d) = sum(T1_nw(:,m,s,d,l)>1.2816)/R;
            size_T1_alrv_10(l,m,s,d) = sum(T1_alrv(:,m,s,d,l)>1.2816)/R;
            size_T1_d_10(l,m,s,d) = sum(T1_d(:,m,s,d,l)>1.2816)/R;
            size_T1_nw_d_10(l,m,s,d) = sum(T1_d_nw(:,m,s,d,l)>1.2816)/R;
            size_T1_alrv_d_10(l,m,s,d) = sum(T1_d_alrv(:,m,s,d,l)>1.2816)/R;
            end
        end
    end
end


str_h1_rowvec = [" "," "," "," "," "," ","h=1"," "," "," "," "," "," "];
str_h4_rowvec = [" "," "," "," "," "," ","h=4"," "," "," "," "," "," "];
str_h12_rowvec = [" "," "," "," "," "," ","h=12"," "," "," "," "," "," "];
str_h24_rowvec = [" "," "," "," "," "," ","h=24"," "," "," "," "," "," "];

str_rho = [" ", "  ","$\rho=0.25"," "," "," ", "$\rho=0.90$", " "," "," ", "$\rho=0.95$", " "," "];

str_mu = ["0.30","0.35","0.40","0.45"];
str_mu_rowvec = ["$\mu_{0}",str_mu,str_mu,str_mu];

str_1 = ["$\beta_{2}=0.0";"$\beta_{2}=0.10";"$\beta_{2}=0.20";"$\beta_{2}=0.30";"$\beta_{2}=0.35";"$\beta_{2}=0.40";"$\beta_{2}=0.45";"$\beta_{2}=0.50";"$\beta_{2}=0.55";"$\beta_{2}=0.60";];

% d versions

T1_d_rho1 = [size_T1_d_10(:,:,1,1);size_T1_d_10(:,:,1,2);size_T1_d_10(:,:,1,3);size_T1_d_10(:,:,1,4)];
T1_d_rho2 = [size_T1_d_10(:,:,2,1);size_T1_d_10(:,:,2,2);size_T1_d_10(:,:,2,3);size_T1_d_10(:,:,2,4)];
T1_d_rho3 = [size_T1_d_10(:,:,3,1);size_T1_d_10(:,:,3,2);size_T1_d_10(:,:,3,3);size_T1_d_10(:,:,3,4)];

T10_d_a = [[str_mu_rowvec];[str_rho];[str_h1_rowvec];[str_1,size_T1_d_10(:,:,1,1),size_T1_d_10(:,:,2,1),size_T1_d_10(:,:,3,1)]];
T10_d_b = [[str_h4_rowvec];[str_1,size_T1_d_10(:,:,1,2),size_T1_d_10(:,:,2,2),size_T1_d_10(:,:,3,2)]];
T10_d_c = [[str_h12_rowvec];[str_1,size_T1_d_10(:,:,1,3),size_T1_d_10(:,:,2,3),size_T1_d_10(:,:,3,3)]];
T10_d_d = [[str_h24_rowvec];[str_1,size_T1_d_10(:,:,1,4),size_T1_d_10(:,:,2,4),size_T1_d_10(:,:,3,4)]];

bigmat_T10_d = [T10_d_a;T10_d_b;T10_d_c;T10_d_d];


%xlswrite('Montec_DGP1a_power.xlsx',bigmat_T10_d,'dgp1a-T10d-power','A1');


 
