function [xtrue] = create_discs(choice,centre1,centre2,r1,r2,n,m,sigma1,sigma2)
% create an image with two circular discs

if choice == 1 % uniform discs
[rr,cc] = meshgrid(1:n);

x1 = sqrt((rr-centre1).^2+(cc-centre1).^2)<r1;
x2 = sqrt((rr-centre2).^2+(cc-centre2).^2)<r2;

% figure, imagesc(xtrue)

else
    if choice == 2 % gaussian discs
x1 = zeros(n,m);

gsize = size(x1);
[R,C] = ndgrid(1:gsize(1), 1:gsize(2));

xc = centre1;
yc = centre1;
exponent = ((R-xc).^2 + (C-yc).^2)./(2*sigma1);
x1      = (exp(-exponent)); 

xc = centre2;
yc = centre2;
exponent = ((R-xc).^2 + (C-yc).^2)./(2*sigma2);
x2      = (exp(-exponent)); 
    end
end

xtrue = x1 + x2;

end



