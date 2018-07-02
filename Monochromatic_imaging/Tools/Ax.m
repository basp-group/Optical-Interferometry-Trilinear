function [Ax] = Ax(x,n,m,c,mask,D)

% This function generates a vector of selected frequency coefficients of x
% considering symmetrisation of the measurement vector

 xobs = uhat_mask(x,n,m,c,mask);

 Ax = D*xobs ;
 
end

