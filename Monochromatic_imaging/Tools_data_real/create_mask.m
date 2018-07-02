function [ M1, M2, M3,count,M_p,M_s_large,M_s_zero,T,mask] = ...
    create_mask(T_uv,n,m,M_b,u_lf,seed,c)

% c : vector position of 0 freq

count = 0;

Ti = T_uv(1:end-1);

[mask_ind,M_p] = generate_idx_mask(n,m, unique(Ti(:)) ) ;


mask = zeros(n*m,1) ;
mask(mask_ind) = 1 ;
% M_p = floor(size(find(mask_ind),1)/2);

M_s_large = reshape(mask,n,m) ;
%figure, imagesc(M_s_large)

M_s = zeros(n,m) ;
M_s(:, 1:m/2) = M_s_large(:, 1:m/2) ;
M_s(1:n/2+1, m/2+1) = M_s_large(1:n/2+1, m/2+1) ;


M_s_conj = zeros(n,m) ;
M_s_conj(:, m/2+2:end) = M_s_large(:, m/2+2:end) ;
M_s_conj(n/2+1:end, m/2+1) = M_s_large(n/2+1:end, m/2+1) ;
M_s_large(n/2+1, m/2+1) = 1 ;

M_s_zero = zeros(n,m) ;
M_s_zero(n/2+1, m/2+1) = 1 ;

M_s = logical(M_s) ;
M_s_conj = logical(M_s_conj) ;
M_s_zero = logical(M_s_zero) ;





% select bispectra -------------------------------------------------------
T = zeros(M_b,3) ;
i0 = 1 ;
% select bispectra from low freq.
M_hf = floor(u_lf*n*m) ;
if M_hf < 2*M_p
[ ~, Tib ] = mask_1D_select_freq( n,m,u_lf, 1/20, seed, c ) ;
Tib = Tib(1:end-1,:) ;
mask_ind_hf = generate_idx_mask(n,m, unique(Tib(:)) ) ;
mask_ind_hf = intersect(mask_ind_hf, mask_ind) ;
mask_hf = zeros(n*m,1) ;
mask_hf(mask_ind_hf) = 1 ;

M_hf = max(size(mask_ind_hf)) ;



while 1
Idx = randperm(M_hf) ;
for j =1:length(Idx)-2
idx = Idx(j:j+2) ;
idx = sort(mask_ind_hf(idx))' ;
if idx(1)~=c && idx(2)~=c && idx(3)~=c
if sum(ismember(idx,T,'rows')) == 0 
if ~(nnz(find(idx(1) == T)) && nnz(find(idx(2) == T)) && nnz(find(idx(3) == T)))
    count = count + 1;
end
T(i0,:) = idx ;
i0 = i0+1 ;
end
end

if i0>floor(0.9*M_b)
    break
end
end
if i0>floor(0.9*M_b) 
    break
end
end


end
% select remaining bispectra
for i = i0:M_b
idx = [c,c,c] ;
while idx(1)==c || idx(2)==c || idx(3)==c
while 1
idx = sort(randperm(M_p*2,3)) ;
idx = mask_ind(idx)' ;
if sum(ismember(idx,T,'rows')) == 0 
    break
end
end
end
if ~(nnz(find(idx(1) == T)) && nnz(find(idx(2) == T)) && nnz(find(idx(3) == T)))
    count = count + 1;
end

T(i,:) = idx ;
end

% build vectors of indices for bispectrum ---------------------------------
M1 = 0*T(:,1) ;
M2 = 0*T(:,2) ;
M3 = 0*T(:,3) ;
for i =1:M_b
M1(i) = sum(mask(1:T(i,1))) ;
M2(i) = sum(mask(1:T(i,2))) ;
M3(i) = sum(mask(1:T(i,3))) ;
end


end

function [mask_ind,M_p] = generate_idx_mask(n,m, mask_ind_inf )

% frequency space sampling: symmetry
mask = zeros(n*m,1) ;
mask(mask_ind_inf) = 1 ;
mask = reshape(mask,n,m) ;
Mask = zeros(n+1,m+1) ;
Mask(1:end-1,1:end-1) = mask ;
Mask = fliplr(flipud(Mask)) ;
mask = Mask(1:end-1,1:end-1) ;

mask = mask(:) ;
mask_ind_sup = find(mask==1) ;

% frequency space sampling
mask_ind = [mask_ind_inf;mask_ind_sup] ;

% number of powerspectra
M_p = length(mask_ind_inf) ;
end

