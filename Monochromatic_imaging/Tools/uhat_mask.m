function [uobs] = uhat_mask(u,n,m,c,mask )

% This function generates a vector of selected frequency coefficients of u

     u = reshape(u,n,m);
     uhat= fft2(u); 
     uhat = fftshift(uhat);
     uhat = uhat(:);
     uobs = uhat(mask);

end

