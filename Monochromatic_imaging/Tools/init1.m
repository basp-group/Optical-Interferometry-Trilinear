function [xi] = init1(y,M_s,param,flag,M_p,M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = param.n*param.m;
c = param.c;

if flag == 1         % initialize with dirty image from power spectrum and actual phases
    xmask = zeros(size(M_s));
 z = sqrt(y(1:M_p)).*exp(ph*sqrt(-1));
 xmask(M_s) = z;
 
   f = flip(xmask);

 xmask(c+1 : N) = conj(f(1: N-c)); 
 xmask(c) = (y(M));

 xmask2 = ifftshift(reshape(xmask,param.n,param.m));
 initx = (param.ni)*ifft2(xmask2);

initx = real(initx);
initx(initx<0) = 0;
   xi = initx;
 
else
    if flag == 2  % initialize with dirty image and random phases
        xmask = zeros(size(M_s));
r = rand(M_p,1);
ph = (2*pi).*r ;
z = sqrt(y(1:M_p)).*exp(ph*sqrt(-1));
xmask(M_s(1:c-1)) = z;
xmask = xmask(:);
xmask = xmask(1:c-1);
 
 f = flip(xmask);

 xmask(c+1 : N) = conj(f(1: N-c)); 
  xmask(c) = 1;

 %figure,imagesc(log(abs(xm))); colorbar
 %figure, imagesc(angle(xm)); colorbar
 
 xmask2 = ifftshift(reshape(xmask,param.n,param.m));
 initx = (param.ni)*ifft2(xmask2);

initx = real(initx);
initx(initx<0) = 0;
   xi = initx;
 
    else
        if flag == 3  % initialize with dirty image -only magnitude
        xmask = zeros(size(M_s));

z = sqrt(y(1:M_p));
 xmask(M_s) = z;

 f = flip(xmask);

 xmask(c+1 : N) = conj(f(1: N-c)); 
 xmask(c) = (y(M));

 xmask2 = ifftshift(reshape(xmask,param.n,param.m));
 initx = (param.ni)*ifft2(xmask2);

initx = real(initx);
initx(initx<0) = 0;
   xi = initx;
        else
    
    
%  f = flip(xmask);
% 
%  xmask(c+1 : N) = conj(f(1: N-c)); 
%  xmask(c) = (y(M+1));
% 
%  xmask2 = ifftshift(reshape(xmask,n,m));
%  initx = ifft2(xmask2);
% 
% initx = real(initx);
% initx(initx<0) = 0;
%    xi = initx;
 
    
    if flag == 4 % initialization using FB imposing reality and positivity
    xmask = zeros(size(M_s));

    z = sqrt(y(1:M_p));
    xmask(M_s) = z;
    f = flip(xmask);
  
   xmask(c+1 : N) = conj(f(1: N-c)); 
   xmask(c) = (y(M));

   xmask2 = ifftshift(reshape(xmask,param.n,param.m));
   initx = (param.ni)*ifft2(xmask2);
   xi = fb_init(y,initx,param,M_p,M_s,M);
    else
        if flag == 5  %random initialisation
            xi =rand(param.n,param.m); 
%             xi =xi(:)./sum(xi(:));
        end        

    end
        end
    end
end
end


