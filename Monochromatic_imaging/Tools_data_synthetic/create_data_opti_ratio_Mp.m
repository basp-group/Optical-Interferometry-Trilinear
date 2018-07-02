function [ y, Mask_large, Mask, Mask_zero, Mask_conj, c, M_p, M_b,T1_s,...
        T2_s,T3_s,count] = create_data_opti_ratio_Mp(n,m,u,seed,xtrue, u_b)

% Create measurements y
% Masks for frequency selection: 
% Mask_large - net mask for u-v plane
% Mask - mask for one half of u-v plane
% Mask_zero - mask to select zero frequency coefficient
% Mask_conj - mask for other half of u-v plane (the conjugate plane)
% T1_s, T2_s, T3_s - selection operators for Fourier coefficients of u1,
%                    u2, and u3, respectively


% zero freq position
xhat= fft2(xtrue); 
xhat = fftshift(xhat);
val = xhat(floor(n/2)+1,floor(m/2)+1);
nxhat = xhat(:);
c = find(nxhat == val);     % index for zero frequency component

% create bispectra
sigmau = 1/4; % standard deviation of the gaussian profile
[T, Ti] = mask_1D_select_freq(n,m,u,sigmau,seed,c);

M_b = floor(u_b*n*m);

% create mask for frequency selection
% and for bispectra in frequency selection

u_lf = 0.02 ; % Ratio of bispectrum points from low frequency region

[Mask, Mask_conj, Mask_zero, Mask_large, M1, M2, M3, M_p,count] = ...
    create_mask_freq_selection_power_bi(Ti,n,m,M_b, u_lf,seed,c) ;

Mask_large = max(Mask_large, Mask_zero) ;
Mask_large = logical(Mask_large);


M = M_p + M_b; 
L = [M1 M2 M3]' ;
L = L(:) ;
[T1,T2,T3,T1_s,T2_s,T3_s] = L1_L2_L3_completeF(M+1,M_p,M_b,L) ;

y = symm_data_2(T1,T2,T3, xtrue,n,m,c,Mask_large) ;


end