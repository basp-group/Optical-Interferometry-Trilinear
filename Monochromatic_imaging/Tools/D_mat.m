function d23 = D_mat(L2_m,L3_m,u2,u3, n,m,c,mask)

% The measurements y correspond to the triple product : T1u1 .* T2u2 .* T3u3
% This function generates the product T2u2.*T3u3, for given u2 and u3


d2 = Ax(u2,n,m,c,mask,L2_m) ;
d3 = Ax(u3,n,m,c,mask,L3_m) ;

d23 = d2.*d3;

end

