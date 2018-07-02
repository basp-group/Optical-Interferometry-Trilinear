function [SNR,u,u1,u2,u3,SNR1,SNR2,SNR3,Time,crit,k] = ...
    solve_blockfb_lambda_hyper_diff_images(xtrue,Y,init,method, L1_m, L2_m, L3_m, ...
    NbIt, NbIt1, NbIt2, NbIt3, param, VM,nu, rw,Psi, Psit,mask_large,l,plots)

%
% Block coordinate forward backward algorithm

% solves for min_{u1,u2,u3} f(u1,u2,u3) + r(u1) + r(u2) + r(u3)

%f(u1,u2,u3) = Data fidelity term: squred l2 norm- smooth function
% r(u_p) : regularisation term for u_p, p = {1,2,3} :
% l2,1 norm + indicator function for positivity and reality
% non- smooth term


% The input argument param contains the following fields:
%
%   - xtrue: Original image (required to compute SNR)
%   - Y : measurements
%   - init : initialisation for u1, u2, u3
%   - method = 0 for only posiitivty, 2 for l1 minimisation + positivity
%   - NbIt : number of global iterations 
%   - NbIt1, NbIt2, NbIt3 : number of forward backward iterations for each
%   variable u1, u2, u3, respectively
%   - Psit: Sparsifying transform 
%   - Psi: Adjoint of Psit
%

% References :
%1) J. Birdi, A. Repetti, and Y. Wiaux. "A regularized tri-linear approach
%  for optical interferometric imaging", submitted to MNRAS, arXiv:
%  1609.00546 [astro-ph.IM].

%2) E. Chouzenoux, J.-C. Pesquet and A. Repetti, ?A block coordinate variable
%  metric forward-backward algorithm,? Accepted in J. Global Optim., 2016.

if rw == 1 % rw = 1 for reweighting
param.Psi = Psi;
param.Psit = Psit;
end

param.nu = nu;
display = 50 ;


U1 = init;
U2 = init;
U3 = init;

u1 =cell(1,l); u2=u1; u3=u1;



for i =1:l
u1{i}=U1(:,:,i);
u2{i}=U2(:,:,i);
u3{i}=U3(:,:,i);
end



u = cell(1,l);
A1 = cell(1,l); A2 = cell(1,l); A3 = cell(1,l);
d23 = cell(1,l); d13 = cell(1,l); d12 = cell(1,l);
Ainv = cell(1,l); Phi1 = cell(1,l); Phi1T = cell(1,l); Au = cell(1,l);
Phi2 = cell(1,l); Phi2T = cell(1,l);
Phi3 = cell(1,l); Phi3T = cell(1,l); SNR = cell(1,l);
SNR1 = cell(1,l) ; SNR2 = cell(1,l) ; SNR3 = cell(1,l) ; 

gamma = 1.9 ;
%  mask_large = param.mask_large ;
n = param.n ;
m = param.m ;
c = param.c ;
C_fft = param.ni ;
if max(size(param.lambda)) == 1
if param.lambda ==1
lambda = ones(6*param.M+1,1) ;
else
Mp = param.Mp ;
lambda = ones(6*param.M+1,1) ;
lambda(6*Mp+1:end-1) = param.lambda ;
end
else
lambda = param.lambda;
end
% Y = Y(:,l);

Y_ = bsxfun(@times,[Y;ones(1,l)],sqrt(lambda)) ;
y_ = mat2cell(Y_,size(Y_,1),[ones(1,l)]);


disp('*************************************************')
disp('PALM algorithm')
 

param_norm1.Psi = param.Psi ;
param_norm1.Psit = param.Psit ; 
param_norm1.max_iter = 200 ;
param_norm1.nu = param.nu ;
param_norm1.n = n ;
param_norm1.m = m ;

k_temp = 1 ;

if method == 1
Psit =@(x) param.Psit(reshape(x,param.n,param.m)) ;
for i =1:l
s{i} = abs(Psit(u1{i}));
s{i} = s{i}(:);
end
S = cell2mat(s);
norms = sqrt(sum(S.^2,2));
norm_l21 = sum(norms(:));
l21term1(k_temp) = nu * norm_l21 ;
else
l21term1(k_temp) = 0 ;
end
l21term2(k_temp) = l21term1(k_temp) ;
l21term3(k_temp) = l21term1(k_temp) ;

 for k = 1:NbIt
    
tic
% Solving for u1 ----------------------------------------------------------
u1_old = u1 ;

if NbIt1 >0
    
for i = 1:l
    
d23{i} = D_mat(L2_m{i},L3_m{i},u2{i},u3{i},n,m,c,mask_large{i});

Phi1{i} =@(x) sqrt(lambda).*d23{i}.*Ax(x,n,m,c,mask_large{i},L1_m{i});
Phi1T{i} =@(y) Atrans(sqrt(lambda).*y,mask_large{i},n,m,C_fft,d23{i},L1_m{i}) ;

if VM == 0  
A1{i} = pow_method_operator(Phi1{i},Phi1T{i},n,m);
elseif VM ==1
A1{i} = abs(d23{i}).^2 ;
end
if k == 1
Au{i} = Phi1{i}(u1{i}) ;
end
end
   
if k ==1
AU = cell2mat(Au);
crit_fid(k_temp) =  0.5*norm(AU(:) - Y_(:), 2)^2; 
crit(k_temp) = crit_fid(k_temp) + l21term1(k_temp) + l21term2(k_temp) + l21term3(k_temp) ;
k_temp = k_temp +1 ;
end

[u1, crit_fid(k_temp),l21term1(k_temp)] = FB_iterations...
(u1, Y_,y_, method, gamma,param.nu, A1, Phi1, Phi1T, NbIt1, param_norm1,VM,l,plots);
l21term2(k_temp) = l21term2(k_temp-1) ;
l21term3(k_temp) = l21term3(k_temp-1) ;
crit(k_temp) = crit_fid(k_temp) + l21term1(k_temp) + l21term2(k_temp) + l21term3(k_temp) ;
k_temp = k_temp+1 ;

end

% -------------------------------------------------------------------------
     
% Solving for u2 ----------------------------------------------------------

u2_old = u2 ;
if NbIt2 >0

for i = 1:l
d13{i} = D_mat(L1_m{i},L3_m{i},u1{i},u3{i},n,m,c,mask_large{i});


Phi2{i} =@(x) sqrt(lambda).*d13{i}.*Ax(x,n,m,c,mask_large{i},L2_m{i});
Phi2T{i} =@(y) Atrans(sqrt(lambda).*y,mask_large{i},n,m,C_fft,d13{i},L2_m{i}) ;
if VM == 0 
A2{i} = pow_method_operator(Phi2{i},Phi2T{i},n,m);
elseif VM ==1
A2{i} = abs(d13{i}).^2 ;
end
end

[u2, crit_fid(k_temp),l21term2(k_temp)] = FB_iterations...
(u2,Y_, y_, method, gamma, param.nu,A2, Phi2, Phi2T, NbIt2, param_norm1,VM,l,plots);
l21term1(k_temp) = l21term1(k_temp-1) ;
l21term3(k_temp) = l21term3(k_temp-1) ;
crit(k_temp) = crit_fid(k_temp) + l21term1(k_temp) + l21term2(k_temp) + l21term3(k_temp) ;

k_temp = k_temp+1 ;
end
% -------------------------------------------------------------------------

% Solving for u3 ----------------------------------------------------------
u3_old = u3 ;

if NbIt3 >0

for i = 1:l
d12{i} = D_mat(L1_m{i},L2_m{i},u1{i},u2{i},n,m,c,mask_large{i});

Phi3{i} =@(x) sqrt(lambda).*d12{i}.*Ax(x,n,m,c,mask_large{i},L3_m{i});
Phi3T{i} =@(y) Atrans(sqrt(lambda).*y,mask_large{i},n,m,C_fft,d12{i},L3_m{i}) ;
if VM == 0 
A3{i} = pow_method_operator(Phi3{i},Phi3T{i},n,m);
elseif VM ==1
A3{i} = abs(d12{i}).^2 ;
end
end

[u3, crit_fid(k_temp),l21term3(k_temp)] = FB_iterations...
(u3,Y_, y_, method, gamma,param.nu, A3, Phi3, Phi3T, NbIt3, param_norm1,VM,l,plots);

l21term1(k_temp) = l21term1(k_temp-1) ;
l21term2(k_temp) = l21term2(k_temp-1) ;
crit(k_temp) = crit_fid(k_temp) + l21term1(k_temp) + l21term2(k_temp) + l21term3(k_temp) ;

k_temp = k_temp+1 ;
end

% -------------------------------------------------------------------------
Time(k) = toc ;   

% updates -----------------------------------------------------------------
nm_u1 = cell(1,l); nm_u2 = cell(1,l); nm_u3 = cell(1,l);

for i = 1:l
nm_u1{i} = norm(u1{i} - u1_old{i})/norm(u1{i}) ;
nm_u2{i} = norm(u2{i} - u2_old{i})/norm(u2{i}) ;
nm_u3{i} = norm(u3{i} - u3_old{i})/norm(u3{i}) ;
SNR1{i}(k) = 20*log(norm(xtrue{i}(:))/norm(xtrue{i}(:)-u1{i}(:)));
SNR2{i}(k) = 20*log(norm(xtrue{i}(:))/norm(xtrue{i}(:)-u2{i}(:)));
SNR3{i}(k) = 20*log(norm(xtrue{i}(:))/norm(xtrue{i}(:)-u3{i}(:)));
end

norm_u1(k) = max(cellfun(@max,nm_u1));
norm_u2(k) = max(cellfun(@max,nm_u2));
norm_u3(k) = max(cellfun(@max,nm_u3));

% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Display -----------------------------------------------------------------
if mod(k, display) == 0
display_results(k,Time,crit,norm_u1,norm_u2,norm_u3) ;
display_figures(crit,u1{1},u2{1},u3{1},norm_u1,norm_u2,norm_u3,SNR1,SNR2,SNR3,n,m) ;
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
 
% Stopping Criterion

d = sort([norm_u1(k),norm_u2(k),norm_u3(k)]);

if method == 0
if k>2 && d(3) <= 1e-03
break
end
else
if method == 1
if d(3) <= 0.005
break
end
end
end
end

disp(['Numberof iterations:',num2str(k)]);

for i=1:l
u{i} = (1/3)*(u1{i} + u2{i} + u3{i});
%     SNR{i} = (1/3)*(SNR1{i}(end) + SNR2{i}(end) + SNR3{i}(end));
SNR{i} = 20*log(norm(xtrue{i}(:))/norm(xtrue{i}(:)-u{i}(:)));
end

end






% *************************************************************************
% *************************************************************************
function [u,crit_fid,crit_pen] = FB_iterations(u,Y,y, method, gamma,nu, A, Phi, Phit, NbIt, param,VM,l,plots)


g = cell(1,l);  dummy = cell(1,l); s = cell(1,l); Au = cell(1,l);
Ainv = cell(1,l); B_dummy = cell(1,l);
B = zeros(l,1); Binv = zeros(l,1);
normg = cell(1,l); norm_it = cell(1,l); 


% param.gamma = param.nu*gamma ;
param.gamma = nu;
for i = 1:l
Ainv{i} = 1./A{i} ;
Au{i} = Phi{i}(u{i}) ;
B(i) = sqrt(gamma*Ainv{i});
Binv(i) = 1./B(i);

if method == 1 
Psit =@(x)  param.Psit(reshape(x,param.n,param.m)) ;
s{i} = abs(Psit(u{i})); 
s{i} = s{i}(:);
end
end

AU = cell2mat(Au);

fidterm(1) =  0.5*norm(AU(:) - Y(:), 2)^2; 

if method == 1 
S = cell2mat(s);
norms = sqrt(sum(S.^2,2));
l21term(1) = nu *sum(norms(:));

e_norm(1) = fidterm(1) + l21term(1);
else
e_norm(1) = fidterm(1) ;
l21term(1) = 0 ;
end


for iter = 1:NbIt
for i = 1:l
% gradient step 
g{i} = Phit{i}(Au{i}-y{i});
g{i} = real(g{i});
normg{i}(iter) = sum(g{i}(:).^2) ;
dummy{i} = u{i}-gamma*Ainv{i}.*g{i};
if method == 1
B_dummy{i} = Binv(i)*dummy{i};
end
end

% proximity step
if method == 1 
dummy = solver_prox_L21real(B_dummy,param,l,B); 
for i = 1:l
dummy{i} = B(i)*dummy{i};
end
end 

% updates  

for i = 1:l
dummy{i} = max(dummy{i},0) ;
norm_it{i}(iter) = sum((u{i}(:)-dummy{i}(:)).^2) ;
Au{i} = Phi{i}(dummy{i}) ;
if method == 1 
Psit =@(x)param.Psit(reshape(x,param.n,param.m)) ;
s{i} = abs(Psit(u{i}));
s{i} = s{i}(:);
end
end

AU = cell2mat(Au);

fidterm(iter+1) =  0.5*norm(AU(:) - Y(:), 2)^2; 

if method == 1 
S = cell2mat(s);
norms = sqrt(sum(S.^2,2));
l21term(iter+1) = param.nu * sum(norms(:));

e_norm(iter+1) = fidterm(iter+1) +l21term(iter+1);
else
e_norm(iter+1) = fidterm(iter+1) ;
l21term(iter+1) = 0 ;  
end


% update
u = dummy;
end

crit_fid = fidterm(end) ;
crit_pen = l21term(end);

if plots
figure(200)
subplot 311, plot(fidterm), xlabel('fidterm')
subplot 312, plot(l21term), xlabel('l21term')
subplot 313, plot(e_norm), xlabel('e_norm')
pause(0.1)
end
end


function display_results(k,Time,crit,norm_u1,norm_u2,norm_u3)
disp(' ')
disp('------------------------------------')
disp(['Iteration : ', num2str(k)])
disp(['Time = ',num2str(sum(Time))])
disp(['Crit = ',num2str(crit(end))])
disp('------------------------------------')
disp('u1 : ')
disp(['|| u_k - u_{k-1} || = ',num2str(norm_u1(end))])
disp('------------------------------------')
disp('u2 : ')
disp(['|| u_k - u_{k-1} || = ',num2str(norm_u2(end))])

disp('------------------------------------')
disp('u3 : ')
disp(['|| u_k - u_{k-1} || = ',num2str(norm_u3(end))])

disp('------------------------------------')
disp('------------------------------------')
end

function display_figures(crit,u1,u2,u3,norm_u1,norm_u2,norm_u3,SNR1,SNR2,SNR3,n,m)
figure(99)
semilogy(crit)
ylabel('crit-with l21')
figure(100)
subplot 231
imagesc(reshape(u1,n,m)), axis image
subplot 232
imagesc(reshape(u2,n,m)), axis image
subplot 233
imagesc(reshape(u3,n,m)), axis image
subplot 234
semilogy(norm_u1), ylabel('u1 : norm evol')
subplot 235
semilogy(norm_u2), ylabel('u2 : norm evol')
subplot 236
semilogy(norm_u3), ylabel('u3 : norm evol')
end
