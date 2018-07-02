function [ut] = Atrans(y,M_s,n,m,C_fft,d,L)

% % % Atrans = F_t * Mf_t * L_t * D_t
% % % with _t --> transpose-conjugate operator
% % %      D = Diag(d)


xt = (conj(d).*y); % size M_p
xt = L'* xt ;

% transposition of mask selecting frequencies
xmask = zeros(size(M_s));
xmask(M_s) = xt;
% % % size(xmask)

xmask2 = ifftshift(reshape(xmask,n,m));
ut = C_fft*ifft2(xmask2);

end

