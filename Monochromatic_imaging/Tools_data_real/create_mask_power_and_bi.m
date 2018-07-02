function [ M_p, M_b,count,M_s_large,M_s_zero,T_bi,mask ] = ...
    create_mask_power_and_bi(n,m,c,seed,u_b,u_lf,T_uv )

% M_p = param.Mp;
M_b = floor(u_b*n*m) ;

% create mask for frequency selection
% and for bispectra in freq. sel.

[ ~, ~, ~,count,M_p,M_s_large,M_s_zero,T_bi,mask] = ...
    create_mask(T_uv,n,m,M_b,u_lf,seed,c) ;



M_s_large = max(M_s_large, M_s_zero);
M_s_large = logical(M_s_large);

end