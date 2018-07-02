
disp([''])
disp(['seed:',num2str(seed)]);
disp(['l21 regularised case, t_max = ',num2str(t_max)]);
disp(['----------------------------'])

init = zeros(n,m,l);

if rec_single_channel
for i =1:l
init(:,:,i) = reshape(sol_l1{i},param.m,param.n);
end
else
for i =1:l
init(:,:,i) = reshape(sol_init{i},param.m,param.n);
end
end

param.nu = mu_l21; 

method = 1;
NbIt1 = t_max ;
NbIt2 = t_max ;
NbIt3 = t_max ;
NbIt = K ;


[SNR,sol,u1sol,u2sol,u3sol,SNR1, SNR2, SNR3,Time, crit,k] = ...
solve_blockfb_lambda_hyper_diff_images(x_bar,Y,init, method, T1_s, T2_s, T3_s, ...
NbIt, NbIt1, NbIt2, NbIt3,param,VM,param.nu,0,param.Psi,param.Psit,Mask_large,l,plots); 
  
SNR_l21{seed} = {cat(1, SNR{:})};
sol_l21{seed} = {cat(1, sol)};
