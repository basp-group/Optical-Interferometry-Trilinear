%% Add noise

noise_choice = 'input_snr' ;

yo = cell(1,l) ;
y = cell(1,l) ;
count = cell(1,l) ; % Number of independent measurements
T1_s= cell(1,l); T2_s= cell(1,l);  T3_s= cell(1,l);
Mask_large = cell(1,l); sigma = cell(1,l);

rw = 0; %Flag for reweighting

disp([''])
disp(['seed:',num2str(seed)]);
disp(['positivity and reality constrained case']);
disp(['----------------------------'])

for ij = 1:l
rng(seed);

[ yo{ij}, M_p{ij}, M_b{ij},T1_s{ij},T2_s{ij},T3_s{ij},count{ij},Mask_large{ij},Mask_zero{ij},c ] = ...
create_data_discrete_grid(n,m,seed,x_bar{ij},u_b,T{ij},Mp_vec{ij},Mp,plots) ;

   
NB=numel(yo{ij});

% Add Gaussian i.i.d. noise
eB=sqrt(1/NB*sum(abs(yo{ij}(1:end)).^2));
switch noise_choice
case 'var'
   param.sigma_noise = 1e-06 ;
case 'input_snr'
    param.sigma_noise= 10^(-iSNR/20)*eB;
end

sigma{ij} = param.sigma_noise;

nBr=(randn(size(yo{ij}))); 
nBi=(randn(size(yo{ij}(6*M_p{ij}+1:end))));


yv = zeros(NB,1);
yv(1:6*M_p{ij}) = yo{ij}(1:6*M_p{ij}) + nBr(1:6*M_p{ij})*sigma{ij};
yv(6*M_p{ij}+1:end) = yo{ij}(6*M_p{ij}+1:end) + (nBr(6*M_p{ij}+1:end) + 1i*nBi)*sigma{ij}/sqrt(2);


y{ij} =yv;   
    

end
 
Y = cell2mat(y); 

 
%% Initialisation

xi = cell(1,I); init = cell(1,I);

for seed_init = 1:I
rng(seed_init)

xi_rand = rand(param.m,param.n);

init{seed_init} = rand(param.m,param.n);
init{seed_init} = init{seed_init}./sum(init{seed_init}(:));


for i =1:l
xi{seed_init}(:,:,i) = xi_rand./sum(xi_rand(:));
end
end



%% Reconstruction with positivity and reality constraints


% Define variables
SNR_old = cell(1,I) ;
sol_old = cell(1,I) ;
u1sol_old = cell(1,I) ; u2sol_old = cell(1,I) ; u3sol_old = cell(1,I) ;
SNR1_old = cell(1,I) ; SNR2_old = cell(1,I) ; SNR3_old = cell(1,I) ;
Time_old = cell(1,I) ; crit_old = cell(1,I) ;
val_old = cell(1,I) ; k_old = cell(1,I) ; DR_old = cell(1,I);
ind = cell(1,l);


val_m = cell(1,I) ; ind = cell(1,l);
sol_init = cell(1,I) ; SNR_init = cell(1,l); 

% 
% 
% % -------------------------------------------------------------------------
% % Algorithm --------------------------------------------------------
% % -------------------------------------------------------------------------

VM = 0 ; 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


NbIt1_old = t_max ; 
NbIt2_old = t_max ;
NbIt3_old = t_max ;
NbIt_old = K; 


%without l1 
method = 0 ;
nu =0;

initialize = 1;

for iter = 1:l
parfor seed_init = 1:I

[SNR_old{seed_init},DR,sol_old{seed_init},...
u1sol_old{seed_init},u2sol_old{seed_init},u3sol_old{seed_init}, ...
SNR1_old{seed_init}, SNR2_old{seed_init}, SNR3_old{seed_init}, ...
Time_old{seed_init}, crit_old{seed_init},k_old{seed_init}] = ...
solve_blockfb_lambda...
(x_bar{iter},y{iter},init{seed_init}, method, T1_s{iter}, T2_s{iter}, T3_s{iter}, ...
NbIt_old, NbIt1_old, NbIt2_old, NbIt3_old,param,VM,nu,0,param.Psi,param.Psit,Mask_large{iter},plots); 

val_old{seed_init} = crit_old{seed_init}(end);

end
 
val_min = min([val_old{:}]);
ind{iter} = find([val_old{:}]==val_min);

sol_init{iter} = sol_old{ind{iter}};
SNR_init{iter} = SNR_old{ind{iter}};
end

SNR_pos{seed} = {cat(1, SNR_init{:})};
x_star{seed} = {cat(1, sol_init{:})};



