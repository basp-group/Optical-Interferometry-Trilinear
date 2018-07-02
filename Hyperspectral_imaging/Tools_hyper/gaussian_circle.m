function [mat] = gaussian_circle(n,m,sigma) 

center = [floor(n/2)+1,floor(m/2)+1];
mat = zeros(n,m);

gsize = size(mat);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));

xc = center(1);
yc = center(2);
exponent = ((R-xc).^2 + (C-yc).^2)./(2*sigma);
mat      = (exp(-exponent)); 
end