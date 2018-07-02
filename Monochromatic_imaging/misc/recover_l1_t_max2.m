
parfor seed = 1:nb_seed
    
disp(['seed:',num2str(seed)]);

disp(['l1 regularised case, t_max = ',num2str(t_max)]);

method = 1 ;
rw = 0;
VM = 0 ;
dwtmode('per');

NbIt1 = t_max ;
NbIt2 = t_max ;
NbIt3 = t_max ;
NbIt = K ;



[SNR_l2{seed},DR_l2{seed},sol_l2{seed},...
u1sol_l2{seed},u2sol_l2{seed},u3sol_l2{seed}, ...
SNR1_l2{seed}, SNR2_l2{seed}, SNR3_l2{seed}, ...
Time_l2{seed}, crit_l2{seed},k_l2{seed}] = ...
solve_blockfb_lambda...
(x_bar,y{seed},int{seed}, method, T1_s{seed}, T2_s{seed}, T3_s{seed}, ...
NbIt, NbIt1, NbIt2, NbIt3, param,VM,mu_l1,rw,param.Psi,param.Psit, Mask_large{seed}); 
 
end