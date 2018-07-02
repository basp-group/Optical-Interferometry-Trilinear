%----------------------------------------------------------------%
% Code for optical interferometric image reconstruction 
% Monochromatic imaging

% uses Block-coordinate forward backward algorithm

% Solves for following constraints:  positivity + reality, l1 with
% different number of sub-iterations, weighted l1, reweighted l1

% Related article :
% J. Birdi, A. Repetti, and Y. Wiaux. "A regularized tri-linear approach
% for optical interferometric imaging", submitted to MNRAS, arXiv:
% 1609.00546 [astro-ph.IM].

%Contact: jb36@hw.ac.uk
%----------------------------------------------------------------%

clc 
clear all
close all

addpath(genpath('misc'));
addpath(genpath('Tools'));
addpath(genpath('Algorithms'));

%% Global parameters

plots = 1; % auxiliary flag for plots
iSNR = 30; % input signal-to-noise ratio
nb_seed = 10; % Number of seeds for different simulations
I = 15; % Number of initialisations
uv_cover = 1; % 1: synthetic uv coverage; 2: realistic uv coverage


%%  Image parameters

% Initial image ---------------------------------------
% LkHa
% Specify the image dimensions
n = 64; % n rows
m = 64; % m columns
N = n*m;  % Size of image vector

x_bar = fitsread('LkHa_ref-img.fits');
x_bar = imresize(x_bar,[n m]);
% -----------------------------------------------------

param.n = n;
param.m = m; 

s = sum(x_bar(:));
x_bar = x_bar/s;       % Normalised true image
if(plots), figure, imagesc(x_bar),colorbar,  axis image, end

%%
if uv_cover == 1
addpath(genpath('Tools_data_synthetic'));
else
addpath('Object1_2016'); 
addpath(genpath('Tools_data_real'));
% addpath('Tools_parallel');
end

%% Create observations

% M_p : Number of power spectrum measurements
% M_b : Number of bispectrum measurements
% M = M_p + M_b : Total number of measurements

u_p = 0.05; % Powerspectrum undersampling ratio, u_p = M_p/N for synthetic coverage
u_b = 0.2; % Bispectrum undersampling ratio, u_b = M_b/N

uf = 2*u_p;  % Fourier undersampling ratio, uf = 2*M_p/N 

seed = 0;

% Create data and undersampling masks

if uv_cover == 1
[ yo, Mask_large, Mask, Mask_zero, Mask_conj, c, M_p, M_b,~,~,~,~ ] = ...
create_data_opti_ratio_Mp(n,m,uf,seed,x_bar,u_b) ;

else
real_cover
end

M = M_p + M_b;  
param.Mp = M_p ;
param.Mb = M_b ;
param.M  = M   ;
param.c = c;
param.ni =1;    % fft coefficient
param.mask_large = Mask_large;  % Under-sampling mask

disp(['Optical interferometric image reconstruction'])
disp(['Monochromatic imaging'])
disp('----------------------------------------')
disp('----------------------------------------')
disp(['number of pixels in the image : ',num2str(n*m)])
disp(['number of Fourier measurements : ',num2str(2*M_p+1)])
disp(['number of measurements : M = ',num2str(M)])
disp(['number of powerspectra, M_p : ',num2str(M_p)])
disp(['number of bispectra, M_b : ',num2str(M_b)])


%% Sparsity operators

% Choose dictionary
% 0 : SARA
% 1 : DB1 == Haar
% 8 : DB8
% 10 : Id
dict = 8;

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
    param.lambda = [ones(size(yo));1];
    param.lambda(6*M_p+1:end-1) = 0 ;
    case 'bi'
    % JUST BI SPECTRUM
    param.lambda = [ones(size(yo));1];
    param.lambda(1:6*M_p) = 0 ;
end

%% Define variables

def_variables

%% Reconstruction with positivity and reality constraints

t_max = 200 ; % Number of FB sub-iterations 
K = 50 ; % Maximum number of global iterations

recover_pos

% x_star  : Solution for 10 different seeds
% SNR_pos : SNR for 10 different seeds

int = x_star; % initialisation for l1/weighted l1 minimisation

%% Reconstruction with l1

t_max = 200; % Number of FB sub-iterations 
t_max1 = t_max;
mu_l1 = 1e-05; % Regularisation parameter for l1

recover_l1_t_max1

% sol_l1 : Solution for 10 different seeds
% SNR_l1 : SNR for 10 different seeds

%% Reconstruction with l1 for diff. t_max


t_max = 400; % Number of FB sub-iterations 
t_max2 = t_max;
mu_l1 = 1e-05; % Regularisation parameter for l1

recover_l1_t_max2

% sol_l2 : Solution for 10 different seeds
% SNR_l2 : SNR for 10 different seeds
 

%% Reconstruction with weighting scheme

t_max = 200; % Number of FB sub-iterations 
mu_rw = 1.5e-05; % Regularisation parameter for weighted l1

recover_l1_weighting_scheme

% sol_r1 : Solution for first weighting iteration for 10 different seeds
% SNR_r1 : SNR for first weighting iteration for 10 different seeds
% sol_r2 : Solution for second weighting iteration for 10 different seeds
% SNR_r2 : SNR for second weighting iteration for 10 different seeds

%%

disp(['SNR_positivity = ',num2str(mean(cellfun(@mean,SNR_pos)))])
disp(['t_max = ',num2str(t_max1), 'SNR_l1 = ',num2str(mean(cellfun(@mean,SNR_l1)))])
disp(['t_max = ',num2str(t_max2), 'SNR_l2 = ',num2str(mean(cellfun(@mean,SNR_l2)))])
disp(['SNR_weightedl1 = ',num2str(mean(cellfun(@mean,SNR_r1)))])
disp(['SNR_reweightedl1 = ',num2str(mean(cellfun(@mean,SNR_r2)))])

