clear all; clc

FONT_SIZE = 20;
FONT_NAME = 'Times';

target = load('target.mat');
target = double(target.data);

prediction = load('prediction.mat');
prediction = double(prediction.data);

%%   >>>>>>>>>>>>>>>>> Let's consider only the smaller domain <<<<<<<<<<<<<<<<<<<<<<
global Soil_profile; global wm; global hm; global w_speed; global G; global rho; global beta; global gamma;
global ne; global nn; global h; global GK; global GM, global GC; global dt; global Meff;
global corner_of_e; global corner_of_b; global corner_of_i; global Ntimestep;
global size_of_i; global size_of_b; global size_of_e; global n_i; global n_b; global n_ge; global time_intg_type;

Soil_profile = 'Material_profile_4';                                % DRM with NO Inclusions
size_of_small_domain = [40,20];
rho = 1500;
h = 1;
dt = 0.003;
Ntimestep = 1.5/0.003;
time_intg_type = 1 ;

wm = size_of_small_domain(1);                                       % width <----> of smaller domain
hm = size_of_small_domain(2);                                       % height of smaller domain

% w_speed = [150;150;150;150];
w_speed = [200;200;200;200];
G = w_speed.^2 .* rho;                                              % shear modulus

% Find elements and nodes
ne = (wm/h)*(hm/h);                                                 % number of elements in the smaller domain
nn = (wm/h+1)*(hm/h+1);                                             % number of nodes in the smaller domain

% Build Global Matrices for Smaller Domain
% if iter == 1
[GK_smaller, GM_smaller, GC_smaller] = FEM_matrices_bigger( );  % build matrices for the SMALLER domain
GK = GK_smaller; GM = GM_smaller; GC = GC_smaller;

Meff_smaller = GM + 0.5 * GC * dt + 0.25 * GK * dt^2;           % Effective matrix to use in the time integration
Meff = Meff_smaller;
% else
% GK = GK_smaller; GM = GM_smaller; GC = GC_smaller; Meff = Meff_smaller;
% end
beta = 0.25; gamma = 0.5;                               % Parameters of Implicit time integration

%%  Find n_i, n_b, n_ge, and n_e again
size_of_e = size_of_small_domain;
corner_of_e = [(wm-size_of_e(1))/2,hm-size_of_e(2)];    % location [x,y] of the bottom-left corner of smaller domain
size_of_b = [size_of_e(1)-10,size_of_e(2)-5];           % [width,height] of Gamma_b
corner_of_b = [corner_of_e(1)+5,corner_of_e(2)+5];      % location [x,y] of the bottom-left corner of Gamma_b
size_of_i = [size_of_b(1)-2,size_of_b(2)-1];            % [width,height] of interior domain (not counting Gamma_b)
corner_of_i = [corner_of_b(1)+1,corner_of_b(2)+1];      % location [x,y] of the bottom-left node of interior domain (not counting Gamma_b)

nodes_small_domain();                                   % defining which nodes are in Gamma_b, Gamma_e, interior domain for the smaller domain
n_b2 = flip(n_b(1:(size_of_b(2)+1)));
n_b2 = [n_b2,n_b((size_of_b(2)+2):end)];

n_ge2 = flip(n_ge(1:(size_of_b(2)+2)));
n_ge2 = [n_ge2,n_ge((size_of_b(2)+3):end)];

%%
for sample = 1:size(target,1)

    target1 = target(sample,:,:);
    target1 = reshape(target1, 126, 501);

    prediction1 = prediction(sample,:,:);
    prediction1 = reshape(prediction1, 126, 501);

    F_eff_target = zeros(nn,Ntimestep+1);
    F_eff_target([n_b2,n_ge2],1:length(target1)) = target1;

    F_eff_predicted = zeros(nn,Ntimestep+1);
    F_eff_predicted([n_b2,n_ge2],1:length(prediction1)) = prediction1;

    %%
    u_history_target = solve_u_DRM2(F_eff_target);

    u_history_predicted = solve_u_DRM2(F_eff_predicted);

    %%
    error_ui(sample,1) = sum(abs(u_history_predicted(n_i,:)-u_history_target(n_i,:)).^2)/sum(abs(u_history_target(n_i,:)).^2)*100;
    % error_ui is the error between targe and predicted wave responses inside the DRM.

    error(sample,1) = sum(abs(F_eff_predicted([n_b,n_ge],:)-F_eff_target([n_b,n_ge],:)).^2)/sum(abs(F_eff_target([n_b,n_ge],:)).^2)*100;

end

filename = 'homo_force_histogram.pdf';

figure;
histogram(error, 10, 'EdgeColor','k', 'FaceColor','y')
xlabel('Error $\% $', 'interpreter','latex','FontSize',FONT_SIZE);
ylabel('Number of samples', 'interpreter','latex','FontSize',FONT_SIZE);
ax = gca;
ax.FontSize = FONT_SIZE;
set(gca,'fontname', FONT_NAME)

set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf, filename, 'Resolution',300)

filename = 'homo_displacement_histogram.pdf';

figure;
histogram(error_ui, 10, 'EdgeColor','k', 'FaceColor','#ffb6c1')
xlabel('Error $\% $', 'interpreter','latex','FontSize',FONT_SIZE);
ylabel('Number of samples', 'interpreter','latex','FontSize',FONT_SIZE);
ax = gca;
ax.FontSize = FONT_SIZE;
set(gca,'fontname', FONT_NAME)

set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf, filename, 'Resolution',300)