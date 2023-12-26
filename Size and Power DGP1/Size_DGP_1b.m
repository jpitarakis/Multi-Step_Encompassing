R=10000;
pi0=0.25;

% dgp 

mu_w = zeros(2,1);
Sigma_w = [1,-0.4;-0.4,0.25];

beta1 = 0.30;
beta2 = 0.0;

mu_vec = [0.30;0.35;0.40;0.45];
    
rho_vec = [0.25;0.90;0.95];
n_vec = [250;500;1000];
h_vec = [1;4;12;24];

[sn,~]=size(n_vec);
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
size_T1_d_nw_10 = nan(sn,sm,sr,sh);
size_T1_d_alrv_10 = nan(sn,sm,sr,sh);

tic

rng = 171966;


for l = 1:sn
n = n_vec(l);

k0 = round(n*pi0);
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
            size_T1_d_nw_10(l,m,s,d) = sum(T1_d_nw(:,m,s,d,l)>1.2816)/R;
            size_T1_d_alrv_10(l,m,s,d) = sum(T1_d_alrv(:,m,s,d,l)>1.2816)/R;
            end
        end
    end
end

mat1_10 = [size_T1_10(:,:,1,1),size_T1_10(:,:,2,1),size_T1_10(:,:,3,1)];
mat2_10 = [size_T1_10(:,:,1,2),size_T1_10(:,:,2,2),size_T1_10(:,:,3,2)];
mat3_10 = [size_T1_10(:,:,1,3),size_T1_10(:,:,2,3),size_T1_10(:,:,3,3)];
mat4_10 = [size_T1_10(:,:,1,4),size_T1_10(:,:,2,4),size_T1_10(:,:,3,4)];

mat1_nw_10 = [size_T1_nw_10(:,:,1,1),size_T1_nw_10(:,:,2,1),size_T1_nw_10(:,:,3,1)];
mat2_nw_10 = [size_T1_nw_10(:,:,1,2),size_T1_nw_10(:,:,2,2),size_T1_nw_10(:,:,3,2)];
mat3_nw_10 = [size_T1_nw_10(:,:,1,3),size_T1_nw_10(:,:,2,3),size_T1_nw_10(:,:,3,3)];
mat4_nw_10 = [size_T1_nw_10(:,:,1,4),size_T1_nw_10(:,:,2,4),size_T1_nw_10(:,:,3,4)];

mat1_alrv_10 = [size_T1_alrv_10(:,:,1,1),size_T1_alrv_10(:,:,2,1),size_T1_alrv_10(:,:,3,1)];
mat2_alrv_10 = [size_T1_alrv_10(:,:,1,2),size_T1_alrv_10(:,:,2,2),size_T1_alrv_10(:,:,3,2)];
mat3_alrv_10 = [size_T1_alrv_10(:,:,1,3),size_T1_alrv_10(:,:,2,3),size_T1_alrv_10(:,:,3,3)];
mat4_alrv_10 = [size_T1_alrv_10(:,:,1,4),size_T1_alrv_10(:,:,2,4),size_T1_alrv_10(:,:,3,4)];


mat1_d_10 = [size_T1_d_10(:,:,1,1),size_T1_d_10(:,:,2,1),size_T1_d_10(:,:,3,1)];
mat2_d_10 = [size_T1_d_10(:,:,1,2),size_T1_d_10(:,:,2,2),size_T1_d_10(:,:,3,2)];
mat3_d_10 = [size_T1_d_10(:,:,1,3),size_T1_d_10(:,:,2,3),size_T1_d_10(:,:,3,3)];
mat4_d_10 = [size_T1_d_10(:,:,1,4),size_T1_d_10(:,:,2,4),size_T1_d_10(:,:,3,4)];

mat1_d_nw_10 = [size_T1_d_nw_10(:,:,1,1),size_T1_d_nw_10(:,:,2,1),size_T1_d_nw_10(:,:,3,1)];
mat2_d_nw_10 = [size_T1_d_nw_10(:,:,1,2),size_T1_d_nw_10(:,:,2,2),size_T1_d_nw_10(:,:,3,2)];
mat3_d_nw_10 = [size_T1_d_nw_10(:,:,1,3),size_T1_d_nw_10(:,:,2,3),size_T1_d_nw_10(:,:,3,3)];
mat4_d_nw_10 = [size_T1_d_nw_10(:,:,1,4),size_T1_d_nw_10(:,:,2,4),size_T1_d_nw_10(:,:,3,4)];

mat1_d_alrv_10 = [size_T1_d_alrv_10(:,:,1,1),size_T1_d_alrv_10(:,:,2,1),size_T1_d_alrv_10(:,:,3,1)];
mat2_d_alrv_10 = [size_T1_d_alrv_10(:,:,1,2),size_T1_d_alrv_10(:,:,2,2),size_T1_d_alrv_10(:,:,3,2)];
mat3_d_alrv_10 = [size_T1_d_alrv_10(:,:,1,3),size_T1_d_alrv_10(:,:,2,3),size_T1_d_alrv_10(:,:,3,3)];
mat4_d_alrv_10 = [size_T1_d_alrv_10(:,:,1,4),size_T1_d_alrv_10(:,:,2,4),size_T1_d_alrv_10(:,:,3,4)];

str_T = ["T=250";"T=500";"T=1000"];
str_T_colvec = [str_T;str_T;str_T;str_T];
T_col = repmat(str_T_colvec,1,1);
str_rho = [" ", "  ","$\rho=0.25"," "," "," ", "$\rho=0.90$", " "," "," ", "$\rho=0.95$", " "," "];
str_mu = ["0.30","0.35","0.40","0.45"];
str_mu_rowvec = ["$\mu_{0}",str_mu,str_mu,str_mu];
str_h1_rowvec = [" "," "," "," "," "," ","h=1"," "," "," "," "," "," "];
str_h4_rowvec = [" "," "," "," "," "," ","h=4"," "," "," "," "," "," "];
str_h12_rowvec = [" "," "," "," "," "," ","h=12"," "," "," "," "," "," "];
str_h24_rowvec = [" "," "," "," "," "," ","h=24"," "," "," "," "," "," "];

heading_1 = [str_rho;str_mu_rowvec;str_h1_rowvec];
heading_4 = [str_h4_rowvec];
heading_12 = [str_h12_rowvec];
heading_24 = [str_h24_rowvec];



% Using known mu0 in variances 

T_10 = [heading_1;str_T,mat1_10;heading_4;str_T,mat2_10;heading_12;str_T,mat3_10;heading_24;str_T,mat4_10];
T_10_nw = [heading_1;str_T,mat1_nw_10;heading_4;str_T,mat2_nw_10;heading_12;str_T,mat3_nw_10;heading_24;str_T,mat4_nw_10];
T_10_alrv = [heading_1;str_T,mat1_alrv_10;heading_4;str_T,mat2_alrv_10;heading_12;str_T,mat3_alrv_10;heading_24;str_T,mat4_alrv_10];

% Using dhat based variances

T_d_10 = [heading_1;str_T,mat1_d_10;heading_4;str_T,mat2_d_10;heading_12;str_T,mat3_d_10;heading_24;str_T,mat4_d_10];
T_d_nw_10 = [heading_1;str_T,mat1_d_nw_10;heading_4;str_T,mat2_d_nw_10;heading_12;str_T,mat3_d_nw_10;heading_24;str_T,mat4_d_nw_10];
T_d_alrv_10 = [heading_1;str_T,mat1_d_alrv_10;heading_4;str_T,mat2_d_alrv_10;heading_12;str_T,mat3_d_alrv_10;heading_24;str_T,mat4_d_alrv_10];

%xlswrite('Montec_DGP1b_size.xlsx',T_10,'dgp1a-chom','A1');
%xlswrite('Montec_DGP1b_size.xlsx',T_10_nw,'dgp1a-nw','A1');
%xlswrite('Montec_DGP1b_size.xlsx',T_10_alrv,'dgp1a-alrv','A1');
%xlswrite('Montec_DGP1b_size.xlsx',T_d_10,'dgp1a-d-chom','A1');
%xlswrite('Montec_DGP1b_size.xlsx',T_d_nw_10,'dgp1a-d-nw','A1');
%xlswrite('Montec_DGP1b_size.xlsx',T_d_alrv_10,'dgp1a-d-alrv','A1');











toc