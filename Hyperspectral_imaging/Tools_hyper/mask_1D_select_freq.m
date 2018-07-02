function [ T, Ti ] = mask_1D_select_freq( n,m,p, sigmau, seed, c )

% mask_1D_select_freq generates a sampling pattern 
%
% The sampling pattern is obtained by sampling frequencies from a 
% bidimensional Gaussian profile in the corresponding Fourier plane. 
% The procedure is repeated until floor(p*n*m) samples are obtained. 
% Note that the originally continuous frequencies are associated with their 
% nearest discrete neighbour,and if a product is sampled twice the result 
% is discarded.
%
% INPUTS
%
% N = n*m dimension of the signal
% p: undersampling ratio
% sigmau: correspond to the stdev of the Gaussian profile
%
% OUTPUT
%
% T table: (floor(p*N) +1) rows (number of sampled frequencies) and 
% floor(p*N) columns - last row to select zero frequency coefficient


rng(seed)

nmeas = n*m*p;

[x,y] = meshgrid(linspace(-1, 1, m), linspace(-1, 1, n)); % Cartesian grid

meas=0; mask=zeros(n,m);
T=zeros(floor(nmeas/2));
Ti=zeros(floor(nmeas/2));
UV = [] ;

while meas<floor(nmeas/2)

% -------------------------------------------------------
% create u
% -------------------------------------------------------
    u=1 ;
    while u == 1 || ui >= c
    while 1
    u=randn(1,2)*sigmau;
    UV = [UV;u] ;
    % find the nearest neighbor on the discrete grid
    d=2;idx_min=0;
    for idx=1:n
        if abs(u(1,1)-x(1,idx))< d
            idx_min=idx; d=abs(u(1,1)-x(1,idx));
        end
    end
    ux=idx_min;

    d=2;idy_min=0;
    for idy=1:n
        if abs(u(1,2)-y(idy,1))< d
            idy_min=idy; d=abs(u(1,2)-y(idy,1));
        end
    end
    uy=idy_min;
    if ux ~=1 && uy ~=1
        break
    end
    end
    
    % find frequency index
    mask(ux,uy)=1; 
    ui = find(mask(:)) ;
    mask=ifftshift(mask);
    u=find(mask(:));
    mask=zeros(n,m);
    end
    

    
    if sum(ismember(T,u))==0
        T(meas+1,:)=u;
        Ti(meas+1,:) = ui ;
        meas=meas+1;
    end
    
end
    
        T(meas+1,:)=1; % open the mask for the F0 frequency
        Ti(meas+1,:) = c ;


% % % figure
% % % plot(UV(:,1), UV(:,2),'.'),axis square

end

