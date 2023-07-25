clear all; clc

FONT_SIZE = 20;
FONT_NAME = 'Times';

category = 'best';
% category = 'q1';
% category = 'median';
% category = 'q3';
% category = 'worst';

error_indices = load('error_indices.mat');

if strcmp(category, 'best')

    sample = error_indices.min(1) + 1

elseif strcmp(category, 'worst')

    sample = error_indices.max(1) + 1

elseif strcmp(category, 'median')

    sample = error_indices.median(1) + 1

elseif strcmp(category, 'q3')

    sample = error_indices.q3(1) + 1

elseif strcmp(category, 'q1')

    sample = error_indices.q1(1) + 1

end

filename = strcat('hetero_', category, '_force_3.pdf');

target = load('target.mat');
target = double(target.data);

prediction = load('prediction.mat');
prediction = double(prediction.data);

target1 = target(sample,:,:);
target1 = reshape(target1, 126, 501);

prediction1 = prediction(sample,:,:);
prediction1 = reshape(prediction1, 126, 501);

% add zeros
F_target = zeros(126, 501);
F_target(:,1:length(target1)) = target1;

F_predicted = zeros(126,501);
F_predicted(:,1:length(prediction1)) = prediction1;

%%
Ntimestep = 500;
load('MyColormap.mat');

%Grid Matrix
xb = [0:60];
xge = [0:64];
yb = [1:Ntimestep+1];
yge = [1:Ntimestep+1];
[xb,yb] = meshgrid(xb,yb);
[xge,yge] = meshgrid(xge,yge);

min_value = min(min(F_target(:, :))); %-6e-9;
max_value = max(max(F_target(:, :))); %6e-9;

if max_value > -min_value

    lower_limit = -max_value;
    upper_limit = max_value;

elseif max_value < -min_value

    lower_limit = min_value;
    upper_limit = -min_value;

end

% figure;
% title('Target Force on $\Gamma_e$','interpreter','latex','FontSize',FONT_SIZE);
% contourf(xge,yge,F_target(62:126,:)',50,'LineColor','none')
% axis off
% % colorbar;
% caxis([lower_limit upper_limit])
% xlabel('$k$ in $P_{e_{kj}}$','interpreter','latex','FontSize',FONT_SIZE)
% ylabel('$j$ in $P_{e_{kj}}$','interpreter','latex','FontSize',FONT_SIZE)

figure_F = figure(1);
subplot(2,2,1)
hold on
title('Target Force on $\Gamma_e$','interpreter','latex','FontSize',FONT_SIZE);
contourf(xge,yge,F_target(62:126,:)',50,'LineColor','none');
% clabel(C, h, 'FontSize', FONT_SIZE)
% colorbar;
caxis([lower_limit upper_limit])
xlabel('$k$ in $P_{e_{kj}}$','interpreter','latex','FontSize',FONT_SIZE)
ylabel('$j$ in $P_{e_{kj}}$','interpreter','latex','FontSize',FONT_SIZE)
plot([0 0],[0 500],'-k')
plot([64 64],[0 500],'-k')
plot([0 64],[500 500],'-k')

ax = gca;
ax.FontSize = FONT_SIZE;

subplot(2,2,2)
hold on
title('Predicted Force on $\Gamma_e$','interpreter','latex','FontSize',FONT_SIZE);
contourf(xge,yge,F_predicted(62:126,:)',50,'LineColor','none')
colorbar;
caxis([lower_limit upper_limit])
xlabel('$k$ in $P_{e_{kj}}$','interpreter','latex','FontSize',FONT_SIZE)
ylabel('$j$ in $P_{e_{kj}}$','interpreter','latex','FontSize',FONT_SIZE)
plot([0 0],[0 500],'-k')
plot([64 64],[0 500],'-k')
plot([0 64],[500 500],'-k')

% min_value = min(min(F_target(1:61, :))); %-6e-9;
% max_value = max(max(F_target(1:61, :))); %6e-9;

ax = gca;
ax.FontSize = FONT_SIZE;

subplot(2,2,3)
hold on
title('Target Force on $\Gamma_b$','interpreter','latex','FontSize',FONT_SIZE);
contourf(xb,yb,F_target(1:61,:)',50,'LineColor','none')
% colorbar;
caxis([lower_limit upper_limit])
% caxis([-0.5 0.5])
xlabel('$k$ in $P_{b_{kj}}$','interpreter','latex','FontSize',FONT_SIZE)
ylabel('$j$ in $P_{b_{kj}}$','interpreter','latex','FontSize',FONT_SIZE)
plot([0 0],[0 500],'-k')
plot([60 60],[0 500],'-k')
plot([0 60],[500 500],'-k')

ax = gca;
ax.FontSize = FONT_SIZE;

subplot(2,2,4)
hold on
title('Predicted Force on $\Gamma_b$','interpreter','latex','FontSize',FONT_SIZE);
contourf(xb,yb,F_predicted(1:61,:)',50,'LineColor','none')
colorbar;
caxis([lower_limit upper_limit])
% caxis([-0.5 0.5])
xlabel('$k$ in $P_{b_{kj}}$','interpreter','latex','FontSize',FONT_SIZE)
ylabel('$j$ in $P_{b_{kj}}$','interpreter','latex','FontSize',FONT_SIZE)
plot([0 0],[0 500],'-k')
plot([60 60],[0 500],'-k')
plot([0 60],[500 500],'-k')

colormap(mymap)

ax = gca;
ax.FontSize = FONT_SIZE;
ax.FontName = FONT_NAME;

set(gcf, 'Position', get(0, 'Screensize'))
exportgraphics(gcf, filename, 'Resolution',300)