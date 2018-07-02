function Dummy = prox(Res,nu)


norms = sqrt(sum(Res.^2,2));
% vec = max(norms-nu, 0)./norms;

vec_num = max(norms-nu,0);
vec = zeros(size(vec_num));
t = logical(norms);

vec(t) = vec_num(t)./norms(t);

Dummy = bsxfun(@times,Res,vec);

% test = sign(Res).*max(abs(Res)-nu, 0);



% C = zeros(size(A));
% t = logical(B);
% C(t) = A(t)./B(t);
