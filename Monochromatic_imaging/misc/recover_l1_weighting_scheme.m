parfor seed = 1:nb_seed
% % weighting l1 : iter 1


disp(['seed:',num2str(seed)]);


rw = 1;
VM = 0 ;
method = 1;
dwtmode('per');

NbIt1 = t_max ;
NbIt2 = t_max ;
NbIt3 = t_max;
NbIt = K ;

sol_n = int{seed};
temp=param.Psit(reshape(sol_n,param.n,param.n));
delta=std(temp(:));

iter = 1;
sigma_n = sigma*sqrt(numel(y{seed})/(numel(x_bar)*1));


delta=max(sigma_n/10,delta);
fprintf('RW iteration: %i\n', iter);
fprintf('delta = %e\n', delta);
init = sol_n;

% Weights
weights=abs(param.Psit1(reshape(sol_n,param.n,param.m)));
weights=delta./(delta+weights);

Psit = @(x) weights.*(param.Psit1(x));
Psi = @(x) param.Psi1(weights.*x);

% Weighted L1 problem
disp(['weighted l1 regularised case, iter = 1']);


[SNR_r1{seed},DR_r1{seed},sol_r1{seed},...
u1sol_r1{seed},u2sol_r1{seed},u3sol_r1{seed}, ...
SNR1_r1{seed}, SNR2_r1{seed}, SNR3_r1{seed}, ...
Time_r1{seed}, crit_r1{seed},k_r1{seed}] = ...
solve_blockfb_lambda...
(x_bar,y{seed},sol_n, method, T1_s{seed}, T2_s{seed}, T3_s{seed}, ...
NbIt, NbIt1, NbIt2, NbIt3,param,VM,mu_rw,rw,Psi,Psit, Mask_large{seed}); 


% Reweighting l1 : iter 2
delta = delta/10;
sol_n = sol_r1{seed};


delta=max(sigma_n/10,delta);

% Weights
weights=abs(param.Psit1(reshape(sol_n,param.n,param.m)));
weights=delta./(delta+weights);

Psit = @(x) weights.*(param.Psit1(x));
Psi = @(x) param.Psi1(weights.*x);

% Reweighted L1 problem
disp(['weighted l1 regularised case, iter = 2']);


[SNR_r2{seed},DR_r2{seed},sol_r2{seed},...
u1sol_r2{seed},u2sol_r2{seed},u3sol_r2{seed}, ...
SNR1_r2{seed}, SNR2_r2{seed}, SNR3_r2{seed}, ...
Time_r2{seed}, crit_r2{seed},k_r2{seed}] = ...
solve_blockfb_lambda...
(x_bar,y{seed},sol_n, method, T1_s{seed}, T2_s{seed}, T3_s{seed}, ...
NbIt, NbIt1, NbIt2, NbIt3,param,VM,mu_rw,rw,Psi,Psit, Mask_large{seed}); 


end

