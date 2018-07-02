function [L1,L2,L3,L1_s,L2_s,L3_s,...
L1_ps,L2_ps,L3_ps, L1_bs,L2_bs,L3_bs,L0] = L1_L2_L3_completeF(M,M_p,Mb,L,Mp_vec,Mp)


L1 = zeros(M,1);
L1(1:Mp) = Mp_vec;
L1(Mp+1:M-1) = L(1:3:end); % bi
L1(end) = M_p+1; % zero

L2 = zeros(M,1);
L2(1:Mp) = Mp_vec;
L2(Mp+1:M-1) = L(2:3:end); % bi
L2(end) = M_p+1; % zero

L3 = zeros(M,1);
L3(1:Mp) = M_p+1; % corresponds to zero frequency
L3(Mp+1:M-1) = L(3:3:end);
L3(end) = M_p+1; % zero


L1_m = zeros(M,2*M_p+1);
L2_m = zeros(M,2*M_p+1);
L3_m = zeros(M,2*M_p+1);
for i = 1:M
    L1_m(i,L1(i))= 1;
    L2_m(i,L2(i))= 1;
    L3_m(i,L3(i))= 1;
end


L1_m =sparse(L1_m);
L2_m =sparse(L2_m);
L3_m =sparse(L3_m);

% powerspectrum
L_p = zeros(Mp, 2*M_p+1) ;

for i = 1:Mp
    L_p(i,L1(i))= 1;
end

L_p_0freq = zeros(size(L_p)) ;
L_p_0freq(:, 1:M_p+1) = L3_m(1:Mp, 1:M_p+1) ;
L_conj = zeros(size(L_p)) ;
L_conj(:,M_p+2:end) = fliplr(L_p(:,1:M_p)) ;


% bispectrums from mask
L1_b = L1_m(Mp+1:end-1, :) ;
L2_b = L2_m(Mp+1:end-1, :) ;
L3_b = L3_m(Mp+1:end-1, :) ;

% zero frequency
L0 = zeros(1,2*M_p+1) ;
L0(M_p+1) = 1 ;



L1_ps = [ L_p ; ...
          L_p ; ...
          L_conj ; ...
          L_conj ; ...
          L_p_0freq ; ...
          L_p_0freq ] ;
L2_ps = [ L_conj ; ...
          L_p_0freq ; ...
          L_p ; ...
          L_p_0freq ; ...
          L_p ; ...
          L_conj ] ;
L3_ps = [ L_p_0freq ; ...
          L_conj ; ...
          L_p_0freq ; ...
          L_p ; ...
          L_conj ; ...
          L_p ] ;

L1_bs = [ L1_b ; ...
          L1_b ; ...
          L2_b ; ...
          L2_b ; ...
          L3_b ; ...
          L3_b ] ;
L2_bs = [ L2_b ; ...
          L3_b ; ...
          L1_b ; ...
          L3_b ; ...
          L1_b ; ...
          L2_b ] ;
L3_bs = [ L3_b ; ...
          L2_b ; ...
          L3_b ; ...
          L1_b ; ...
          L2_b ; ...
          L1_b ] ;

L1 = [ L1_ps ; L1_bs] ;
L2 = [ L2_ps ; L2_bs ] ;
L3 = [ L3_ps ; L3_bs ] ;
L1_s = [ L1 ; L0 ] ;
L2_s = [ L2 ; L0 ] ;
L3_s = [ L3 ; L0 ] ;


 

end

