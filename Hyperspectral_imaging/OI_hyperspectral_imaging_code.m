%----------------------------------------------------------------%
% Code for optical interferometric image reconstruction 
% Hyperspectral imaging

% uses Block-coordinate forward backward algorithm

% Related article :
% J. Birdi, A. Repetti, and Y. Wiaux. "A regularized tri-linear approach
% for optical interferometric imaging", submitted to MNRAS, arXiv:
% 1609.00546 [astro-ph.IM].

%Contact: jb36@hw.ac.uk

%----------------------------------------------------------------%


clc
clear all
close all

addpath('Object1_2016');
addpath('Tools_hyper');
addpath('Algorithms');
addpath('misc');

%% Global parameters

plots = 0; %auxiliary flag for plots
iSNR = 30; %input signal-to-noise ratio
nb_seed = 10; % Number of seeds for different simulations
I = 15; % Number of initialisations

%%  Image parameters

%Initial image ----------------------------------------
% 1: LkHa; 
% 2: discs 
s = 1;

% Specify the image dimensions
n = 64; % n rows
m = 64; % m columns
N = n*m;  % Size of image vector


switch s
case 1
x_bar = fitsread('LkHa_ref-img.fits');
x_bar = imresize(x_bar,[n m]);

case 2
choice = 1; % 1: uniform disc: 2 : Gaussian
centre1 =15;  r1 = 5;
centre2 = 40;  r2 = 10;
sigma1 = 5; sigma2 = 15;
x_bar = create_discs(choice,centre1,centre2,r1,r2,n,m,sigma1,sigma2);
x_bar = x_bar.*10;
end


param.n = n;
param.m = m; 

x_bar=(x_bar+abs(x_bar))./2;

s = sum(x_bar(:));
x_bar = x_bar/s;       % Normalised true image

if(plots), figure, imagesc(x_bar),colorbar,  axis image, end

%%
disp(['Optical interferometric image reconstruction'])
disp(['Hyperspectral imaging']);
disp('----------------------------------------')
disp('----------------------------------------')

%% Load u-v points

l = 2; % Number of spectral channels

load_uv

 %% Discrete coverage
discretize_coverage

%% Generate images and measurement for each channel

% M_p : Number of power spectrum measurements
% M_b : Number of bispectrum measurements
% M = M_p + M_b : Total number of measurements

seed = 0;

s_len = l-1;

gen_images_diff_lambda;

x_bar = cell(1,l);

u_b = 0.1; % Bispectrum undersampling ratio, u_b = M_b/N

for i = 1 : l
x_bar{i} = xf(:,:,i)./sum(sum(xf(:,:,i)));

a =  x_bar{i};
XF1(:,i) = a(:);

[ yo{i}, ~, ~,~,~,~,~,~,~,c ] = ...
create_data_discrete_grid(n,m,seed,x_bar{i},u_b,T{i},Mp_vec{i},Mp,plots) ;


% [yo{i}, ~, M_b{i},~,~,~,~,~,~,c ] = ...
%     create_data_discrete_grid(n,m,seed,xtrue,u_b,T_uv,Mp_vec,Mp,plots )

end


%%
center_x = n/2;
center_y = m/2;

offset = 10;

rec = zeros(n,m);
rec (center_y-offset:center_y+offset,center_x-offset:center_x+offset) = 1;
rec = rec(:);

B = XF1(logical(rec),:);

%%

param.Mb = floor(u_b*n*m) ; % Number of bispectrum measurements in each channel
param.M  = param.Mp + param.Mb  ; % Total number of measurements in each channel
param.c = c;
param.eta = 1.1;
param.weights=1;
param.ni =1;    % fft coefficient
param.xmax = s; % incident flux
param.zero =1;

%% Sparsity operators

% Choose dictionary
% 0 : SARA
% 1 : DB1 == Haar
% 8 : DB8
% 10 : Id
dict = 10;

choose_dictionary
 
%%
% lambda parameter to scale measurement vector yo
% full  : both bispectrum and power spectrum measurements
% power : only power spectrum measurements
% bi    : only bispectrum measurements


lambda_choice = 'full' ;

switch lambda_choice
    case 'full'
    param.lambda = 1;
    case 'power'
    % JUST POWER SPECTRUM
    param.lambda = [ones(size(yo{1}));1];
    param.lambda(6*M_p+1:end-1) = 0 ;
    case 'bi'
    % JUST BI SPECTRUM
    param.lambda = [ones(size(yo{1}));1];
    param.lambda(1:6*M_p) = 0 ;
end


%% Define variables

def_variables_hyper_spectral


%% Reconstruction

rec_single_channel = 1; % Reconstruction with l1 - Single band reconstruction

for seed = 1:nb_seed
% Reconstruction with positivity and reality constraints

t_max = 200 ; % Number of FB sub-iterations 
K = 50 ; % Maximum number of global iterations

recover_pos_hyper

% x_star  : Solution for 10 different seeds
% SNR_pos : SNR for 10 different seeds

%----------------------------------------------------------------%
%----------------------------------------------------------------%

if rec_single_channel  % Reconstruction with l1 - Single band reconstruction
    
t_max = 200; % Number of FB sub-iterations 
t_maxl1 = t_max;
K = 50 ; % Maximum number of global iterations
mu_l1 = 2e-05; % Regularisation parameter for l1

recover_single_channel

% sol_l11 : Solution for 10 different seeds
% SNR_l11 : SNR for 10 different seeds

end

%----------------------------------------------------------------%
%----------------------------------------------------------------%

% Reconstruction with l2,1 - Hyperspectral cube

t_max = 200 ; % Number of FB sub-iterations 
t_maxl21 = t_max ;
K = 50 ; % Maximum number of global iterations
mu_l21 = 2e-05;  % Regularisation parameter

recover_hyperspectral_cube

% sol_l21 : Solution for 10 different seeds
% SNR_l21 : SNR for 10 different seeds

end



%%
% To analyze results

output_hyper

