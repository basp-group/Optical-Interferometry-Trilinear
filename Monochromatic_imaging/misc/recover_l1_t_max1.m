
parfor seed = 1:nb_seed
% %with l1
method = 1 ;
rw = 0;
VM = 0 ;
disp(['seed:',num2str(seed)]);

disp(['l1 regularised case, t_max = ',num2str(t_max)]);


dwtmode('per');

NbIt1 = t_max ;
NbIt2 = t_max ;
NbIt3 = t_max ;
NbIt = K ;



[SNR_l1{seed},DR_l1{seed},sol_l1{seed},...
u1sol_l1{seed},u2sol_l1{seed},u3sol_l1{seed}, ...
SNR1_l1{seed}, SNR2_l1{seed}, SNR3_l1{seed}, ...
Time_l1{seed}, crit_l1{seed},k_l1{seed}] = ...
solve_blockfb_lambda...
(x_bar,y{seed},int{seed}, method, T1_s{seed}, T2_s{seed}, T3_s{seed}, ...
NbIt, NbIt1, NbIt2, NbIt3,param,VM,mu_l1,rw,param.Psi,param.Psit, Mask_large{seed}); 

end