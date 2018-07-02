function [yo,L1_s,L2_s,L3_s,Mask_large,M,Mp,Mb] = generate_mask_and_data....
                                (n,m,c,seed,u_b,T,l,xtrue)

% Generates sampling masks and measurement vector yo
% Output:
% M : total number of measurements
% Mp : Number of power spectrum measurements
% Mb : Number of bispectrum measurements
% Mask_large : sampling mask in 2D Fourier plane

count = cell(1,l) ; % Number of independent measurements
M_s_large = cell(1,l); 
M_s_zero = cell(1,l);
M_p= cell(1,l); M_b= cell(1,l);
T_bi = cell(1,l);
mask = cell(1,l);
u_lf = 0.02 ; % Proportion of bispectrum measurements from low frequency region



for ij = 1:l

%     dwtmode('per');
[ M_p{ij}, M_b{ij},count{ij},M_s_large{ij},M_s_zero{ij},T_bi{ij},mask{ij}]= ...
   create_mask_power_and_bi(n,m,c,seed,u_b,u_lf,T{ij} );


end   
 
% %% Concatenate masks and measurements for all spectral channels
% 
Mask_large = M_s_large{1};
mask_re = reshape(mask{1},n,m);
T_b = T_bi{1};
T_pp = T{1};

if l>1
for i = 2:l

Mask_large = max(Mask_large,M_s_large{i});

mask_re = max(mask_re,reshape(mask{i},n,m));

T_b = [T_b;T_bi{i}];
T_pp = [T_pp;T_p{i}];

end
end

mask = mask_re(:);

Mb = size(T_b,1);
M_p_short = floor(size(find(Mask_large),1)/2);
Mp = size(T_pp,1)-l;
M = Mp + Mb;



[yo,L1_s,L2_s,L3_s] = create_data_and_bi_mask(T_b,M,Mp,Mb,Mask_large,mask,....
                        xtrue,n,m,c,T_pp,M_p_short);

end

