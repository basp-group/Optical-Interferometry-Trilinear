function val = pow_method_operator(A,At,n,m)
%Computes the maximum eigen value of the compund 
%operator AtA
%   

% % % % A =@(x) D * M * fft2(x) ;
% % %  A=@(x) D*Ax(x,n,m,c,mask,L);
% % % % At =@(x) ifft2( Mt * D' *x ) ;
% % %  At = @(x) Atrans(x,mask,n,m,c,C_fft,D,L) ;
x=randn(n,m);
x=x/norm(x(:));

p = 1 + 10^(-6) ;
pnew = 1 ;

n = 1 ;

epsilon = 10^(-8) ;

nmax = 200;

cond = abs( pnew-p ) / pnew ;

% Iterations

while ( cond >= epsilon && n < nmax)
     xnew=At(A(x));
    p=pnew;
    pnew=norm(xnew) / norm(x) ;
    
    cond = abs(  pnew-p ) / pnew ;
    
    x = xnew;
    n = n+1 ;
% % %     disp(['p: ', num2str(p)])
    
end
val = p ;
% x=randn(im_size);
% x=x/norm(x(:));
% init_val=1;
% 
% for k=1:max_iter
%     x=At(A(x));
%     val=norm(x(:));
%     rel_var=abs(val-init_val)/init_val;
%     if (verbose > 0)
%         fprintf('Iter = %i, norm = %e \n',k,val);
%     end
%     if (rel_var < tol)
%         break;
%     end
%     init_val=val;
%     x=x/val;
%     
% end


end

