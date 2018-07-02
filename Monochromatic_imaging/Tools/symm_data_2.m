function ys = symm_data_2(L1,L2,L3, x,n,m,c,mask)

% This function produces symmetrised version of the measurements.

d1 = Ax(x,n,m,c,mask,L1) ;
d2 = Ax(x,n,m,c,mask,L2) ;
d3 = Ax(x,n,m,c,mask,L3) ;


ys = d1.*d2.*d3 ;


end