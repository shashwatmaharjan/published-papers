% 2D DRM SH wave - randomizer - written by Bruno Guidio - Feb/2023 -
% modified May/2023 - modified May/31/2023

clear all; close all;
tic
global wm; global hm; global h; global rho; global G; global freq; global ne; global nn;
global time_intg_type; global dt; global beta; global gamma;
global Ntimestep;global GK; global GM; global Soil_profile; global GC; global Meff;
global w_speed; global s_loc; global um;
global corner_of_b; global size_of_b; global corner_of_i; global size_of_i;
global corner_of_e; global size_of_e;
global n_i; global n_b; global n_e; global n_ge; global n_e_out;

LASTN = maxNumCompThreads(1);

% DEFINE PARAMETERS
Soil_profile = 'Material_profile_6_bigger';

% Set size of bigger domain
wm = 400;                                               % width of bigger domain
hm = 120;                                               % height of bigger domain
size_of_bigger_domain = [wm,hm];

h = 1;                                                  % element size

size_of_small_domain = [40,20];                         % [width,height] of smaller domain

x_sensors = [5:1:(size_of_small_domain(1)-5)];          % location of sensors in the smaller domain' surface (meters)
s_loc = x_sensors +((size_of_small_domain(1)/h+1)...
    *(size_of_small_domain(2)/h)+1);            % location of sensors in the smaller domain (nodes)

rho = 1500;                                             % mass density of the whole computational domain

% w_speed_bigger = [1800;1500;300;250;200;150];           % wave speed for each layer - from bottom to top layer
w_speed_bigger = [1800;1500;300;300;250;200];
w_speed_smaller = [w_speed_bigger(5);w_speed_bigger(5);...
    w_speed_bigger(6);w_speed_bigger(6)]; % wave speed inside the smaller domain !! ALWAYS DOUBLE CHECK - it changes depending on the size of smaller domain

Ntimestep = 500;                                       % number of timesteps

time_intg_type = 1;                                     % 1 means the implicit time integration (this code only supports implicit t.i)
dt = 3e-3 ;% 0.001;
beta = 0.25; gamma = 0.5;                               % Parameters of Implicit time integration
N_iteration = 1;

for iter = 1:N_iteration

    iter
    wm = size_of_bigger_domain(1);                          % width of bigger domain
    hm = size_of_bigger_domain(2);                          % height of bigger domain
    w_speed = w_speed_bigger;
    G = w_speed.^2 .* rho;                                  % shear modulus

    % Find elements and nodes in the bigger domain
    ne = (wm/h)*(hm/h);                                     % number of elements in the bigger domain
    nn = (wm/h+1)*(hm/h+1);                                 % number of nodes in the bigger domain
    if iter == 1
        [GK_bigger, GM_bigger, GC_bigger] = FEM_matrices_bigger_06_28( );                  % create GK, GC, and GM for the bigger domain
        GK = GK_bigger; GM = GM_bigger; GC = GC_bigger;
    else
        GK = GK_bigger; GM = GM_bigger; GC = GC_bigger;
    end

    %% 
    N_force_min = 1; N_force_max = 2;
    N_force = randi([N_force_min,N_force_max],1,1);
    GF_m = zeros(nn,Ntimestep+1);

    %% 
    for iter_force = 1:N_force

        %%  RANDOMIZER - AMPLITUDE, FREQUENCY, DELAY OF FORCE, LOCATION OF FORCE
        A_min = 50; A_max = 150;
        A = randi([A_min,A_max],1,1);                                                   % Random - Ricker's amplitude

        freq_min = 1; freq_max = 5;                                    % Random - frequency
        freq = round(freq_min + (freq_max-freq_min)*rand(1,1),1);

        delay_min = -0.05; delay_max = 0.0;% 0.5;        % BG:May-31 = max from 0.5 to 0.0, min from 0 to -0.05
        delay = round(delay_min + (delay_max-delay_min)*rand(1,1),3);                   % Random - delay of force

        hm_point_min = 20; hm_point_max = 110;          % BG:May-31 = min from 10 to 20
        hm_point = randi([hm_point_min,hm_point_max],1,1);

        wm_point_min = 50; wm_point_max = 200;          % BG:May-31 = min from 20 to 50
        if hm_point > 90
            wm_point_max = 160;
        end
        wm_point = randi([wm_point_min,wm_point_max],1,1);

        point_load_loc = [wm_point,hm_point];                                           % Random - location

        for i = 1:Ntimestep+1
            t = (i-1)*dt;
            F_m_i(i,1) = ricker2(freq,A,t,delay);
        end

        GF_m((wm/h+1)*(point_load_loc(2)/h)+1+(point_load_loc(1)/h),:) = F_m_i' +...
            GF_m((wm/h+1)*(point_load_loc(2)/h)+1+(point_load_loc(1)/h),:);

    end

    %% 
    size_of_e = size_of_small_domain;                       % [width,height] of smaller domain
    corner_of_e = [(wm-size_of_e(1))/2,hm-size_of_e(2)];    % location [x,y] of the bottom-left corner of smaller domain

    corner_of_b = [corner_of_e(1)+5,corner_of_e(2)+5];      % location [x,y] of the bottom-left corner of Gamma_b
    size_of_b = [size_of_e(1)-10,size_of_e(2)-5];           % [width,height] of Gamma_b

    corner_of_i = [corner_of_b(1)+1,corner_of_b(2)+1];      % location [x,y] of the bottom-left node of interior domain (not counting Gamma_b)
    size_of_i = [size_of_b(1)-2,size_of_b(2)-1];            % [width,height] of interior domain (not counting Gamma_b)

    F_eff_smaller = solve_F_eff_bigger (GF_m);              % Force vector for the smaller domain
    
    %%  >>>>>>>>>>>>>>>>>>>>>>>>>>>>> 2nd PART <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % >>>>>>>>>>>>>>>>> Let's consider only the smaller domain <<<<<<<<<<<<<<<<<<<<<<

    Soil_profile = 'Material_profile_4';                    % DRM with NO Inclusions

    wm = size_of_e(1);                                      % width <----> of smaller domain
    hm = size_of_e(2);                                      % height of smaller domain

    w_speed = 0; w_speed = w_speed_smaller;
    G = 0; G = w_speed.^2 .* rho;                           % shear modulus

    % Find elements and nodes
    ne = (wm/h)*(hm/h);                                     % number of elements in the smaller domain
    nn = (wm/h+1)*(hm/h+1);                                 % number of nodes in the smaller domain

    % Build Global Matrices for Smaller Domain
    if iter == 1
        [GK_smaller, GM_smaller, GC_smaller] = FEM_matrices_bigger_06_28( );  % build matrices for the SMALLER domain
        GK = GK_smaller; GM = GM_smaller; GC = GC_smaller;
        Meff_smaller = GM + 0.5 * GC * dt + 0.25 * GK * dt^2;           % Effective matrix to use in the time integration
        Meff = Meff_smaller;
    else
        GK = GK_smaller; GM = GM_smaller; GC = GC_smaller; Meff = Meff_smaller;
    end
    %%  Find n_i, n_b, n_ge, and n_e again
    corner_of_e = [(wm-size_of_e(1))/2,hm-size_of_e(2)];    % location [x,y] of the bottom-left corner of smaller domain
    corner_of_b = [corner_of_e(1)+5,corner_of_e(2)+5];      % location [x,y] of the bottom-left corner of Gamma_b
    corner_of_i = [corner_of_b(1)+1,corner_of_b(2)+1];      % location [x,y] of the bottom-left node of interior domain (not counting Gamma_b)

    nodes_small_domain();                                   % defining which nodes are in Gamma_b, Gamma_e, interior domain for the smaller domain

    %% 
    um = solve_u_DRM(F_eff_smaller);                        % measurement data at the sensors

    % [u_curt_history_small] = solve_u_DRM2 (F_eff_smaller);

    %%  forces to plot & error
    [F_eff_ge,F_eff_b]  = force_eff_on_gamma (F_eff_smaller);
    F_eff = [F_eff_b;F_eff_ge];                             % Effective forces Gamma_b and Gamma_e

    %%  Save um and F_eff for each iteration
    um_history{iter,1} = um;
    F_eff_history{iter,1} = F_eff;

end

save('F_history.mat', 'F_eff_history', '-v7.3')
save('u_history.mat', 'um_history', '-v7.3')

toc