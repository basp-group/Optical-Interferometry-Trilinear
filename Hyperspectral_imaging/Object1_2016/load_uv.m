% Load u-v points

u_pow = load('ucoord');
v_pow = load('vcoord');

u1_bi = load('u1coord');
v1_bi = load('v1coord');

u2_bi = load('u2coord');
v2_bi = load('v2coord');

lambda = load('eff_wave');

% separation vector b between two telescopes

u_p = u_pow.ucoord(:);     
v_p = v_pow.vcoord(:);

u1_b = u1_bi.u1coord(:);
v1_b = v1_bi.v1coord(:);
u2_b = u2_bi.u2coord(:);
v2_b = v2_bi.v2coord(:);

w = lambda.eff_wave(:);
w = w.*10^6;


choose_wavelength = 1;

if choose_wavelength % wavelength range to be decided separately
w = size(l,1);
w(1) = 1.95 ;
w(l) = 1.97;
end

step_size = (w(l)-w(1))/l;

for i = 2:l-1
    w(i) = w(i-1)+step_size;
end




up_lambda = cell(1,l);  vp_lambda = cell(1,l);
u1_lambda = cell(1,l);  v1_lambda = cell(1,l);
u2_lambda = cell(1,l);  v2_lambda = cell(1,l);
u3_lambda = cell(1,l);  v3_lambda = cell(1,l);


%% 
for i = 1:l
    up_lambda{i} = [u_p./w(i)];
    vp_lambda{i} = [v_p./w(i)];
    u1_lambda{i} = [u1_b./w(i)];
    v1_lambda{i} = [v1_b./w(i)];
    u2_lambda{i} = [u2_b./w(i)];
    v2_lambda{i} = [v2_b./w(i)];  
    
    u3_lambda{i} = -u1_lambda{i}-u2_lambda{i};
    v3_lambda{i} = -v1_lambda{i}-v2_lambda{i};
end

up_l = up_lambda{1};
vp_l = vp_lambda{1};
u1_l = u1_lambda{1};
v1_l = v1_lambda{1};
u2_l = u2_lambda{1};
v2_l = v2_lambda{1};
u3_l = u3_lambda{1};
v3_l = v3_lambda{1};


if l > 1
    for i =2:l
    up_l=[up_l;up_lambda{i}];
    vp_l=[vp_l;vp_lambda{i}];
    u1_l=[u1_l;u1_lambda{i}];
    v1_l=[v1_l;v1_lambda{i}];
    u2_l=[u2_l;u2_lambda{i}];
    v2_l=[v2_l;v2_lambda{i}];
    u3_l=[u3_l;u3_lambda{i}];
    v3_l=[v3_l;v3_lambda{i}];
    end
end

  
%%


step_size = (w(l)-w(1))/l;

for i = 2:l-1
    w(i) = w(i-1)+step_size;
end
