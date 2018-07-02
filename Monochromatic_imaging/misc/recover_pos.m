%% Add noise 

NB = numel(yo);

noise_choice = 'iSNR' ;

% Add Gaussian i.i.d. noise
eB=sqrt(1/NB*sum(abs(yo(1:end)).^2));
switch noise_choice
case 'var'
   param.sigma_noise = 1e-06 ;
case 'iSNR'
    param.sigma_noise= 10^(-iSNR/20)*eB;
end

% noise bound based on a chi-squared model
param.epsilon=sqrt(6)*sqrt(2*NB + 4*sqrt(NB))*param.sigma_noise ;
% tolerance on the noise bound
param.epsilon_up= sqrt(6)*(sqrt(param.sigma_noise^2*(2*NB) + 5*param.sigma_noise^2*sqrt(2*(2*NB))));%expectation + 2.1 std
param.epsilon_low= sqrt(6)*(sqrt(param.sigma_noise^2*(2*NB)  -5*param.sigma_noise^2*sqrt(2*(2*NB))));

sigma = param.sigma_noise;

parfor seed = 1:nb_seed

disp(['seed:',num2str(seed)]);
disp(['positivity and reality constrained case']);


if uv_cover == 1
[ yo, Mask_large{seed}, Mask, Mask_zero, Mask_conj, c, M_p, M_b,T1_s{seed},T2_s{seed},T3_s{seed},count_ind ] = ...
create_data_opti_ratio_Mp(n,m,uf,seed,x_bar,u_b) ;
count{seed} = count_ind ;

else
if uv_cover == 2
[yo,T1_s{seed},T2_s{seed},T3_s{seed},Mask_large{seed},M,M_p,M_b] = generate_mask_and_data....
(n,m,param.c,seed,u_b,T_p,l,x_bar);
                            
end
end


rw = 0; %Flag for reweighting

nBr=(randn(size(yo))); 
nBi=(randn(size(yo(6*M_p+1:end))));

yv = zeros(NB,1);
yv(1:6*M_p) = yo(1:6*M_p) + nBr(1:6*M_p)*sigma;
yv(6*M_p+1:end) = yo(6*M_p+1:end) + (nBr(6*M_p+1:end) + 1i*nBi)*sigma/sqrt(2);

y{seed} =yv;    
    
 
%% Initialisation

xi = cell(1,I);

for seed_init = 1:I
rng(seed_init)
xi{seed_init} = rand(param.n,param.m); 
end


%% Reconstruction with positivity and reality constraints

% Define variables
SNR_old = cell(1,I) ;
sol_old = cell(1,I) ;
u1sol_old = cell(1,I) ; u2sol_old = cell(1,I) ; u3sol_old = cell(1,I) ;
SNR1_old = cell(1,I) ; SNR2_old = cell(1,I) ; SNR3_old = cell(1,I) ;
Time_old = cell(1,I) ; crit_old = cell(1,I) ;
val_old = cell(1,I) ; k_old = cell(1,I) ; DR_old = cell(1,I);

val_m = cell(1,I) ;
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


for seed_init =  1:I
in =  xi{seed_init}./sum(xi{seed_init}(:)); 

disp(['number of initialisation:',num2str(seed_init)]);

%without l1
method = 0;
nu =0;

[SNR_old{seed_init},DR_old{seed_init},sol_old{seed_init},...
u1sol_old{seed_init},u2sol_old{seed_init},u3sol_old{seed_init}, ...
SNR1_old{seed_init}, SNR2_old{seed_init}, SNR3_old{seed_init}, ...
Time_old{seed_init}, crit_old{seed_init},k_old{seed_init}] = ...
solve_blockfb_lambda...
(x_bar,y{seed},in, method, T1_s{seed}, T2_s{seed}, T3_s{seed}, ...
NbIt_old, NbIt1_old, NbIt2_old, NbIt3_old, param,VM,nu,rw,param.Psi,param.Psit, Mask_large{seed}); 

val_old{seed_init} = crit_old{seed_init}(end);

end
%  
val_min = min([val_old{:}]);
ind = find([val_old{:}]==val_min);
  

x_star{seed} = sol_old{ind};
SNR_pos{seed} = SNR_old{ind};
DR_pos{seed} = DR_old{ind};
Time_pos{seed} = Time_old{ind};
crit_pos{seed} = val_old{ind};
k_pos{seed} = k_old{ind};

end


