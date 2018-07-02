function [SNR,DR,u,u1,u2,u3,SNR1,SNR2,SNR3,Time,crit,k] = ...
solve_blockfb_lambda(xtrue,y,init,method, T1, T2, T3, ...
NbIt, NbIt1, NbIt2, NbIt3, param, VM,nu, rw,Psi, Psit,mask_large )


%
% Block coordinate forward backward algorithm

% solves for min_{u1,u2,u3} f(u1,u2,u3) + r(u1) + r(u2) + r(u3)

%f(u1,u2,u3) = Data fidelity term: squred l2 norm- smooth function
% r(u_p) : regularisation term for u_p, p = {1,2,3} - l1 norm
% non- smooth term

% The input argument param contains the following fields:
%
%   - xtrue: Original image (required to compute SNR)
%   - y : measurements
%   - init : initialisation for u1, u2, u3
%   - method = 0 for only posiitivty, 2 for l1 minimisation + positivity
%   - NbIt : number of global iterations 
%   - NbIt1, NbIt2, NbIt3 : number of forward backward iterations for each
%     variable u1, u2, u3, respectively
%   - Psit: Sparsifying transform 
%   - Psi: Adjoint of Psit
%

% References :
%1) J. Birdi, A. Repetti, and Y. Wiaux. "A regularized tri-linear approach
%  for optical interferometric imaging", submitted to MNRAS, arXiv:
%  1609.00546 [astro-ph.IM].

%2) E. Chouzenoux, J.-C. Pesquet and A. Repetti, "A block coordinate variable
%  metric forward-backward algorithm", Accepted in J. Global Optim., 2016.


if rw == 1 % rw = 1 for reweighting
param.Psi = Psi;
param.Psit = Psit;
end

param.nu = nu;
display = 50 ;

u1 = init(:);
u2 = init(:);
u3 = init(:);

gamma = 1.9 ;

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
y_ = sqrt(lambda).*[y;1] ;


disp('*************************************************')

param_norm1.Psi = param.Psi ;
param_norm1.Psit = param.Psit ; 
param_norm1.max_iter = 200 ;
param_norm1.n = n ;
param_norm1.m = m ;

k_temp = 1 ;

if method == 1
Psit =@(x) param.Psit(reshape(x,param.n,param.m)) ;
s = abs(Psit(u1));
l1term1(k_temp) = nu * sum(s(:)) ;
else
l1term1(k_temp) = 0 ;
end
l1term2(k_temp) = l1term1(k_temp) ;
l1term3(k_temp) = l1term1(k_temp) ;

for k = 1:NbIt

tic
% Solving for u1 ----------------------------------------------------------
u1_old = u1 ;
if NbIt1 >0
d23 = D_mat(T2,T3,u2,u3,n,m,c,mask_large );

% Phi large
Phi1 =@(x) sqrt(lambda).*d23.*Ax(x,n,m,c,mask_large,T1);
Phi1T =@(y) Atrans(sqrt(lambda).*y,mask_large,n,m,C_fft,d23,T1) ;
if k == 1
crit_fid(k_temp) = 0.5*sum(abs(Phi1(u1)-y_).^2) ;
crit(k_temp) = crit_fid(k_temp) + l1term1(k_temp) + l1term2(k_temp) + l1term3(k_temp) ;
k_temp = k_temp +1 ;
end
if VM == 0 

A1 = pow_method_operator(Phi1,Phi1T,n,m);
elseif VM ==1
A1 = abs(d23).^2 ;
end
[u1, crit_fid(k_temp),l1term1(k_temp)] = FB_iterations...
(u1, y_, method, gamma,param.nu, A1, Phi1, Phi1T, NbIt1, param_norm1,VM);
l1term2(k_temp) = l1term2(k_temp-1) ;
l1term3(k_temp) = l1term3(k_temp-1) ;
crit(k_temp) = crit_fid(k_temp) + l1term1(k_temp) + l1term2(k_temp) + l1term3(k_temp) ;
k_temp = k_temp+1 ;
end
% -------------------------------------------------------------------------

% Solving for u2 ----------------------------------------------------------
u2_old = u2 ;
if NbIt2 >0
d13 = D_mat(T1,T3,u1,u3,n,m,c,mask_large);
% Phi large
Phi2 =@(x) sqrt(lambda).*d13.*Ax(x,n,m,c,mask_large,T2);
Phi2T =@(y) Atrans(sqrt(lambda).*y,mask_large,n,m,C_fft,d13,T2) ;
if VM == 0
A2 = pow_method_operator(Phi2,Phi2T,n,m);
elseif VM ==1
A1 = abs(d23).^2 ;
end
[u2, crit_fid(k_temp),l1term2(k_temp)] = FB_iterations...
(u2, y_, method, gamma,param.nu, A2, Phi2, Phi2T, NbIt2, param_norm1,VM);
l1term1(k_temp) = l1term1(k_temp-1) ;
l1term3(k_temp) = l1term3(k_temp-1) ;
crit(k_temp) = crit_fid(k_temp) + l1term1(k_temp) + l1term2(k_temp) + l1term3(k_temp) ;
k_temp = k_temp+1 ;
end
% -------------------------------------------------------------------------

% Solving for u3 ----------------------------------------------------------
u3_old = u3 ;
if NbIt3 >0
d12 = D_mat(T1,T2,u1,u2,n,m,c,mask_large );
% Phi large
Phi3 =@(x) sqrt(lambda).*d12.*Ax(x,n,m,c,mask_large,T3);
Phi3T =@(y) Atrans(sqrt(lambda).*y,mask_large,n,m,C_fft,d12,T3) ;
if VM == 0
A3 = pow_method_operator(Phi3,Phi3T,n,m);
elseif VM ==1
A1 = abs(d23).^2 ;
end
[u3, crit_fid(k_temp),l1term3(k_temp)] = FB_iterations...
(u3, y_, method, gamma,param.nu, A3, Phi3, Phi3T, NbIt3, param_norm1,VM);
l1term1(k_temp) = l1term1(k_temp-1) ;
l1term2(k_temp) = l1term2(k_temp-1) ;
crit(k_temp) = crit_fid(k_temp) + l1term1(k_temp) + l1term2(k_temp) + l1term3(k_temp) ;
k_temp = k_temp+1 ;
end
% -------------------------------------------------------------------------
Time(k) = toc ;   

% updates -----------------------------------------------------------------
norm_u1(k) = norm(u1 - u1_old)/norm(u1) ;
norm_u2(k) = norm(u2 - u2_old)/norm(u2) ;
norm_u3(k) = norm(u3 - u3_old)/norm(u3) ;
SNR1(k) = 20*log(norm(xtrue(:))/norm(xtrue(:)-u1));
SNR2(k) = 20*log(norm(xtrue(:))/norm(xtrue(:)-u2));
SNR3(k) = 20*log(norm(xtrue(:))/norm(xtrue(:)-u3));
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Display -----------------------------------------------------------------
if mod(k, display) == 0
display_results(k,Time,crit,norm_u1,SNR1,norm_u2,SNR2,norm_u3,SNR3) ;
display_figures(crit,u1,u2,u3,norm_u1,SNR1,norm_u2,SNR2,norm_u3,SNR3,n,m) ;
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Stopping Criterion

d = sort([norm_u1(k),norm_u2(k),norm_u3(k)]);

if method == 0
if d(3) <= 1e-02
break
end
else
if method == 1
if d(3) <= 1e-02 
break
end
end
end

end
disp(['Number of iterations:',num2str(k)]);
u = (1/3)*(u1+u2+u3);
SNR = calculate_snr(xtrue,reshape(u,param.n,param.m));

%Dynamic Range calculation
d23 = D_mat(T2,T3,u,u,n,m,c,mask_large );
Phi1 =@(x) sqrt(lambda).*d23.*Ax(x,n,m,c,mask_large,T1);
Au = Phi1(u);

Phi1T =@(y) Atrans(sqrt(lambda).*y,mask_large,n,m,C_fft,d23,T1) ;
A1 = pow_method_operator(Phi1,Phi1T,n,m);
DR_num = sqrt(n*m)*A1*max(u(:));
DR_den = Phi1T(Au-y_);
DR_den = norm(DR_den(:));
DR = DR_num/DR_den;

end

% *************************************************************************
% *************************************************************************
function [u,crit_fid,crit_pen] = FB_iterations(u, y, method, gamma,nu, A, Phi, Phit, NbIt, param,VM)

Ainv = 1./A ;
param.gamma = nu*gamma*Ainv ;


Au = Phi(u) ;
fidterm(1) = 0.5*sum(abs(Au-y).^2) ;
if method == 1 
Psit =@(x) param.Psit(reshape(x,param.n,param.m)) ;
s = abs(Psit(u));
l1term(1) = nu * sum(s(:)) ;
e_norm(1) = fidterm(1) + l1term(1) ;
else
e_norm(1) = fidterm(1) ;
l1term(1) = 0 ;
end



for iter = 1:NbIt

% gradient step 
g = Phit(Au-y);
g = real(g(:));
dummy = u-gamma*Ainv.*g;

% proximity step
if method == 1 
dummy = solver_prox_L1orig(reshape(dummy,param.n,param.m),param);
dummy = dummy(:) ;
end
dummy = max(real(dummy),0) ;

% updates     
Au = Phi(dummy) ; 

fidterm(iter+1) = 0.5*sum(abs(Au-y).^2) ;
if method == 1
s = abs(Psit(dummy));
l1term(iter+1) = nu * sum(s(:)) ;
e_norm(iter+1) = fidterm(iter+1) + l1term(iter+1) ;
else
e_norm(iter+1) = fidterm(iter+1) ;
l1term(iter+1) = 0 ;
end


% update
u = dummy;
end

crit_fid = fidterm(end) ;
crit_pen = l1term(end) ;

end


function display_results(k,Time,crit,norm_u1,SNR1,norm_u2,SNR2,norm_u3,SNR3)
disp(' ')
disp('------------------------------------')
disp(['Iteration : ', num2str(k)])
disp(['Time = ',num2str(sum(Time))])
disp(['Crit = ',num2str(crit(end))])
disp('------------------------------------')
disp('u1 : ')
disp(['|| u_k - u_{k-1} || = ',num2str(norm_u1(end))])
disp(['SNR = ',num2str(SNR1(end))])
disp('------------------------------------')
disp('u2 : ')
disp(['|| u_k - u_{k-1} || = ',num2str(norm_u2(end))])
disp(['SNR = ',num2str(SNR2(end))])
disp('------------------------------------')
disp('u3 : ')
disp(['|| u_k - u_{k-1} || = ',num2str(norm_u3(end))])
disp(['SNR = ',num2str(SNR3(end))])
disp('------------------------------------')
disp('------------------------------------')
end

function display_figures(crit,u1,u2,u3,norm_u1,SNR1,norm_u2,SNR2,norm_u3,SNR3,n,m)
figure(99)
semilogy(crit)
ylabel('crit')
figure(100)
subplot 331
imagesc(reshape(u1,n,m)), axis image
subplot 332
imagesc(reshape(u2,n,m)), axis image
subplot 333
imagesc(reshape(u3,n,m)), axis image
subplot 334
semilogy(norm_u1), ylabel('u1 : norm evol')
subplot 335
semilogy(norm_u2), ylabel('u2 : norm evol')
subplot 336
semilogy(norm_u3), ylabel('u3 : norm evol')
subplot 337
plot(SNR1), ylabel('u1 : SNR')
subplot 338
plot(SNR2), ylabel('u2 : SNR')
subplot 339
plot(SNR3), ylabel('u3 : SNR')
pause(0.1)
end
