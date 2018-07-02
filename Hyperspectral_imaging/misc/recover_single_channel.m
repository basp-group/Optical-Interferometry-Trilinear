
NbIt1 = t_max ;
NbIt2 = t_max ;
NbIt3 = t_max ;
NbIt = K ;


%with l1 : each band separately
method = 1 ;
rw = 0;

disp([''])
disp(['seed:',num2str(seed)]);
disp(['l1 regularised case, t_max = ',num2str(t_max)]);
disp(['----------------------------'])

parfor iter = 1:l

dwtmode('per');

[SNR_l1{iter},DR,sol_l1{iter},...
u1sol_l1{iter},u2sol_l1{iter},u3sol_l1{iter}, ...
SNR1_l1{iter}, SNR2_l1{iter}, SNR3_l1{iter}, ...
Time_l1{iter}, crit_l1{iter},k_l1{iter}] = ...
solve_blockfb_lambda...
(x_bar{iter},y{iter},reshape(sol_init{iter},param.m,param.n), method, T1_s{iter}, T2_s{iter}, T3_s{iter}, ...
NbIt, NbIt1, NbIt2, NbIt3,param,VM,mu_l1,0,param.Psi,param.Psit,Mask_large{iter}); 

end


SNR_l11{seed} = {cat(1, SNR_l1{:})};
sol_l11{seed} = {cat(1, sol_l1)};

