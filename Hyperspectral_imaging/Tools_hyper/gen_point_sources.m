function [im] = gen_point_sources(n,m,k)

im1 = zeros(n-20,m-20);
ind = ceil(rand(1,k)*(n-20)*(m-20));
im1(ind) = 1;

im = padarray(im1,[10,10]);

end
