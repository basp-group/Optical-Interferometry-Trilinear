function [y,L1_s,L2_s,L3_s] = create_data_and_bi_mask(T_b,M,M_p,M_b,M_s_large,mask,....
                        xtrue,n,m,c,T_p,M_p_short)


% build vectors of indices for bispectrum ---------------------------------
M1 = 0*T_b(:,1) ;
M2 = 0*T_b(:,2) ;
M3 = 0*T_b(:,3) ;
for i =1:M_b
M1(i) = sum(mask(1:T_b(i,1))) ;
M2(i) = sum(mask(1:T_b(i,2))) ;
M3(i) = sum(mask(1:T_b(i,3))) ;
end

L = [M1 M2 M3]' ;
L = L(:) ;

[L1,L2,L3,L1_s,L2_s,L3_s] = L1_L2_L3_completeF_real(M_p,M_b,L,M_p_short) ;

y = symm_data_2(L1,L2,L3,xtrue,n,m,c,M_s_large) ;

end