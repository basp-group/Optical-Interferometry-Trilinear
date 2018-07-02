
nlevel = 4 ;  % Wavelet decomposition depth
dwtmode('per');
 
switch dict
case 0
disp('SARA basis')
[param.Psi1, param.Psit1] = SARA_basis(x_bar, nlevel) ;
case 1
%DB2 wavelet basis
disp('Haar basis')
[C,S]=wavedec2(x_bar,nlevel,'db1'); 
ncoef=length(C);
param.Psit1 = @(x) wavedec2(x, nlevel,'db1'); 
param.Psi1 = @(x) waverec2(x,S,'db1');
case 8 
%DB8 wavelet basis
disp('DB8 basis')
[C,S]=wavedec2(x_bar,nlevel,'db8'); 
ncoef=length(C);
param.Psit1 = @(x) wavedec2(x, nlevel,'db8'); 
param.Psi1 = @(x) waverec2(x,S,'db8');
case 10
%image
disp('sparsity on the image')
param.Psit1 = @(x) x; 
param.Psi1 = @(x) x;
end


param.Psi = param.Psi1;
param.Psit = param.Psit1;
