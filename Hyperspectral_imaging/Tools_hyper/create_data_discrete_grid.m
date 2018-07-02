function [ y, M_p, M_b,L1_s,L2_s,L3_s,count,M_s_large,M_s_zero,c ] = ...
    create_data_discrete_grid(n,m,seed,xtrue,u_b,T_uv,Mp_vec,Mp,plots )

% zero frequency position

xhat= fft2(xtrue); 
xhat = fftshift(xhat);

val = xhat(floor(n/2)+1,floor(m/2)+1);
nxhat = xhat(:);
c = find(nxhat == val);     % index for zero frequency component
param.c = c;

u_lf = 0.02 ; 

M_b = floor(u_b*n*m) ;

% create mask for frequency selection
% and for bispectra in frequency selection
[ M1, M2, M3,count,M_p,M_s_large,M_s_zero] = ...
create_mask(T_uv,n,m,M_b, u_lf,seed,c,plots) ;



M_s_large = max(M_s_large, M_s_zero);
M_s_large = logical(M_s_large);

M = Mp + M_b;
L = [M1 M2 M3]' ;
L = L(:) ;
[L1,L2,L3,L1_s,L2_s,L3_s] = L1_L2_L3_completeF(M+1,M_p,M_b,L,Mp_vec,Mp) ;

y = symm_data_2(L1,L2,L3, xtrue,n,m,c,M_s_large) ;


end


