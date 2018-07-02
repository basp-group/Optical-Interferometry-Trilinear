
% im : True image at reference frequency
% s : number of channels


disp '/*----------------------------------------------------------------*/'
disp ' Simulation of multi-wavelength optical interferometric data'
disp ' Physical model: power law spectra'
disp '/*----------------------------------------------------------------*/'

untitled4
untitled3

im = x_bar;
s = s_len+1;
%%

[Nix, Niy]=size(im);%Image dimensions

rng(seed);

%%
G = fspecial('gaussian',[10 10],2);

z4=z4./max(z4(:));

z4 = imresize(z4,[64,64]);


beta_im=imfilter(im,G,'same');
beta_im=beta_im./max(beta_im(:));

beta2= 3*beta_im+z4;

% beta2=beta2./max(beta2(:));
beta2 = 0.5*beta2;

if plots
figure(2),
imagesc(beta2), colorbar, axis image, colormap 'jet',
end

%%

xf=zeros(Nix,Niy,length(l));

%%
XF1=zeros(Nix*Niy,l);

for t=1:l

%     xf1(:,:,t)=im .* (f(t)./f0).^(-beta2 + (beta1.*log(f(t)./f0)));
    xf1(:,:,t)=im .* (w(1)./w(t)).^(-beta2);

    a = xf1(:,:,t);
    XF1(:,t) = a(:);
 
end

[m,n,c] = size(xf1);
p = m*n;
[M,N]=size(XF1);

%%
[u,s,v]=svd(XF1, 'econ');
s = diag(s);
r1 = rank(XF1)

%%
XF = XF1;
xf = xf1;
r = r1;


%%


center_x = Nix/2;
center_y = Niy/2;

offset = 20;

rec = zeros(Nix,Niy);
rec (center_y-offset:center_y+offset,center_x-offset:center_x+offset) = 1;
rec = rec(:);

B = XF(logical(rec),:);

if(plots)
for i = 1 :20: size(B,1)
    figure(11);
    plot(B(i,:)),hold on;
end
end
% axis([1 64 0 4]);

%%
if plots
figure,imagesc(log10(XF+1e-4)),colormap jet
figure(5),
subplot 121,imagesc(log10(xf(:,:,1)+1e-4)), colormap 'jet',
subplot 122,imagesc(log10(xf(:,:,end)+1e-4)), colormap 'jet',
%%
figure,
imagesc(log10(im+1e-4)), colormap 'jet',
end
%%
% center_x = Nix/2;
% center_y = Niy/2;
% 
% figure(5),
% imagesc(indmap), axis square, colorbar, colormap 'jet',hold on
% box off;set(gcf, 'Position', [0 400 800 800]);
% % rectangle('Position',[center_x-64 center_x-64 128 128],'LineWidth',2,'EdgeColor','r');
% rectangle('Position',[center_x-10 center_x-10 20 20],'LineWidth',2,'EdgeColor','r');
% 
% axis([0 256 0 256]);
%%
N = Nix * Niy;
Nx = Nix;
Ny = Niy;
c = length(f);




