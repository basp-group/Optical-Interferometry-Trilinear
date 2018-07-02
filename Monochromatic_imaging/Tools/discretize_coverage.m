
nx = param.n; %size of image
 
b_max = cell(1,l); 

pos = cell(1,l); posD = cell(1,l);
ud = cell(1,l); vd = cell(1,l);
mask_ = cell(1,l); mask_e = cell(1,l);
v_rep = cell(1,l); T_p = cell(1,l);
mask_t = cell(1,l);

for i = 1:l
b_max{i} = 2*ceil(2*max(sqrt(up_lambda{i}.^2+vp_lambda{i}.^2)));
end

bmax = max(cellfun(@max,b_max));
du = floor(bmax/nx);
 
    
for i = 1:l
 
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
mask_t{i} = zeros(nx,nx);
mask_t{i}(ud{i}(ii),vd{i}(ii))=1 ;
T_p{i}(ii) = find(mask_t{i}(:));
end
mask_e{i} = reshape(flipud(mask_{i}(:)),nx,nx) ;
T_p{i} = sort(T_p{i}');


%% 
    
T= cell(1,l);
    
for i = 1:l
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

Mp = size(unique(T{i}(:)),1);

param.Mp = Mp;

Mask_l = max(mask_mod,reshape(flipud(mask_mod(:)),nx,nx));

M_s_zero = zeros(param.n,param.m) ;
M_s_zero(param.n/2+1, param.m/2+1) = 1 ;
M_s_zero = logical(M_s_zero) ;



M_p_orig = size(T_p{i}(:),1);
T_p{i}(end+1,:) = c;
T{i}(size(T{i},1)+1,:) = c;

end
end
