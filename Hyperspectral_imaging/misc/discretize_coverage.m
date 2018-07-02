
nx = param.n; %size of the image
     
b_max = cell(1,l); 
pos = cell(1,l); posD = cell(1,l);
ud = cell(1,l); vd = cell(1,l);
mask_ = cell(1,l); mask_e = cell(1,l);
rep_rows = cell(1,l); mask_rep = cell(1,l);
rep_count = cell(1,l); idx_rep = cell(1,l);
num_rep = cell(1,l); count_rep = cell(1,l);

% zero frequency position

xhat= fft2(x_bar); 
xhat = fftshift(xhat);

val = xhat(floor(n/2)+1,floor(m/2)+1);
nxhat = xhat(:);
c = find(nxhat == val);     % index for zero frequency component
param.c = c;


for i = 1:(l)
b_max{i} = 2*ceil(2*max(sqrt(up_lambda{i}.^2+vp_lambda{i}.^2)));
end

bmax = max(cellfun(@max,b_max));
du = floor(bmax/nx);
 
    
for i = 1:(l)

for j=1:length(up_lambda{i})
    pos{i}(j)=up_lambda{i}(j);
    posD{i}(j)=floor(pos{i}(j)/du)*du;
    ud{i}(j)=posD{i}(j);
    pos{i}(j)=vp_lambda{i}(j);
    posD{i}(j)=floor(pos{i}(j)/du)*du;
    vd{i}(j)=posD{i}(j);
end
ud{i} = ud{i}/bmax+ (1/2);
vd{i} = vd{i}/bmax+ (1/2);
ud{i} = ceil(ud{i}*nx)+1 ; vd{i} = ceil(vd{i}*nx) ;

% discrete coverage is given by ud-vd
mask_{i} = zeros(nx,nx) ;
for ii = 1:length(ud{i})
mask_{i}(ud{i}(ii),vd{i}(ii))=1 ;
end
mask_e{i} = reshape(flipud(mask_{i}(:)),nx,nx) ;

% To find repeating points in the discrete coverage
t = [ud{i}(:),vd{i}(:)];
B = sortrows(t);
S = [1;any(diff(B),2)];
[L,S] = regexp(sprintf('%i',S'),'1(0)+','start','end');
rep_rows{i} = B(S,:); % Repeated Rows.
rep_count{i} = (S-L+1)'; % How often each repeated row appears
count_rep{i} = sum(rep_count{i}); % Total count of repeated elements
num_rep{i} = size(rep_rows{i},1); % Number of repeating elements


idx_rep{i} = zeros(length(rep_rows{i}),1);

for ii = 1:length(rep_rows{i})
mask_rep = zeros(nx,nx);
mask_rep(rep_rows{i}(ii,1),rep_rows{i}(ii,2)) = 1;
idx_rep{i}(ii) = find(mask_rep(:));
end
 
end
    

%% 

T= cell(1,l); yo= cell(1,l);
M_p= cell(1,l); M_b= cell(1,l);
L1_s= cell(1,l); L2_s= cell(1,l);  L3_s= cell(1,l);
count= cell(1,l); Mask_large= cell(1,l);
Mask_zero= cell(1,l); Mp_vec = cell(1,l);
    
for i = 1:(l)
mask_mod = mask_{i};
mask_mod(1,:) = 0;
mask_mod(:,1)=0;

% To make sure conjugate pairs are not selected

mask_odd = mask_mod(2:end,2:end);
mm = mask_odd(:);
c_odd = floor(length(mm)/2)+1;
m_first = mm(1:c_odd-1);
m_sec = mm(c_odd+1:end);
m_conj = flipud(m_sec);
m_tab = [m_first,m_conj];
dummy =ones(size(m_tab));
f = ismember(m_tab,dummy,'rows');
rep_ind = find(f); % indices for selected conjugate pairs

if rep_ind
    m_conj(rep_ind) = 0;
end

mask_inv = [m_first;0;flipud(m_conj)];
mask_inv = reshape(mask_inv,n-1,m-1);

mask = zeros(n,m) ;

mask(2:end,2:end) = mask_inv; 


T{i} = find(mask(:));

ff = find(ismember(T{i},idx_rep{i}));

Mp_vec{i} = ones(size(ud{i},2),1);

Mp_mat = zeros(size(ud{i},2)-(count_rep{i}-num_rep{i}));
Mp_mat(1,:) = 1:size(ud{i},2)-((count_rep{i}-num_rep{i}));

for ii = 1:size(ff,1)
    Mp_mat(1:rep_count{i}(ii),ff(ii)) = ff(ii);
end

Mp_mat = Mp_mat(:);
Mp_vec{i} = Mp_mat(find(Mp_mat));


Mp = size(Mp_vec{i},1);


param.Mp = Mp;

Mask_l = max(mask_mod,reshape(flipud(mask_mod(:)),nx,nx));

M_s_zero = zeros(param.n,param.m) ;
M_s_zero(param.n/2+1, param.m/2+1) = 1 ;
M_s_zero = logical(M_s_zero) ;


T{i}(size(T{i},1)+1,:) = c;
end