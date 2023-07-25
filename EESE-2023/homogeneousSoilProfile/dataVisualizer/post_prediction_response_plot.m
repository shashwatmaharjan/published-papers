% Run the response_solver.m first

global h; global hm; global wm
i = 320;
j = i+10;
k = j+10;

xm = [5:h:35]; ym = [5:h:hm];

x = [0:h:wm]; y = [0:h:hm];
[x,y] = meshgrid(x,y);

%%
u_curt = u_history_target(:,i);
W = zeros(length(u_curt),1);
W_prev = u_curt;
W = W_prev;
W1 = reshape(W,[(wm/h+1),(hm/h+1)])';

u_curt = u_history_target(:,j);
W = zeros(length(u_curt),1);
W_prev = u_curt;
W = W_prev;
W2 = reshape(W,[(wm/h+1),(hm/h+1)])';

u_curt = u_history_target(:,k);
W = zeros(length(u_curt),1);
W_prev = u_curt;
W = W_prev;
W3 = reshape(W,[(wm/h+1),(hm/h+1)])';

%
u_curt = u_history_predicted(:,i);
W = zeros(length(u_curt),1);
W_prev = u_curt;
W = W_prev;
W4 = reshape(W,[(wm/h+1),(hm/h+1)])';

u_curt = u_history_predicted(:,j);
W = zeros(length(u_curt),1);
W_prev = u_curt;
W = W_prev;
W5 = reshape(W,[(wm/h+1),(hm/h+1)])';

u_curt = u_history_predicted(:,k);
W = zeros(length(u_curt),1);
W_prev = u_curt;
W = W_prev;
W6 = reshape(W,[(wm/h+1),(hm/h+1)])';

%%
min_value_1 = min(min(W1(6:21,6:36)));
max_value_1 = max(max(W1(6:21,6:36))); %6e-9;

min_value_2 = min(min(W2(6:21,6:36)));
max_value_2 = max(max(W2(6:21,6:36))); %6e-9;

min_value_3 = min(min(W3(6:21,6:36)));
max_value_3 = max(max(W3(6:21,6:36))); %6e-9;

min_value = min([min_value_1, min_value_2, min_value_3]);
max_value = min([max_value_1, max_value_2, max_value_3]);

if max_value > -min_value

    lower_limit = -max_value;
    upper_limit = max_value;

elseif max_value < -min_value

    lower_limit = min_value;
    upper_limit = -min_value;

end

filename = strcat('homo_', category, '_target_response.pdf');

fig1 = figure ('Position',[200 10 850 400]) ;  %680

subplot(2,3,1);
hold on;
title('0.32 s','interpreter','latex','FontSize',FONT_SIZE);
contourf(xm,ym,W1(6:21,6:36),50,'LineColor','none')
ylabel('$y$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
xlabel('$x$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
% colorbar;
caxis([min_value max_value])

ax = gca;
ax.FontSize = FONT_SIZE;
set(gca,'fontname', FONT_NAME)

subplot(2,3,2);
hold on;
title('0.33 s','interpreter','latex','FontSize',FONT_SIZE);
contourf(xm,ym,W2(6:21,6:36),50,'LineColor','none')
ylabel('$y$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
xlabel('$x$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
% colorbar;
caxis([min_value max_value])

ax = gca;
ax.FontSize = FONT_SIZE;
set(gca,'fontname', FONT_NAME)

subplot(2,3,3);
hold on;
title('0.34 s','interpreter','latex','FontSize',FONT_SIZE);
contourf(xm,ym,W3(6:21,6:36),50,'LineColor','none')
ylabel('$y$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
xlabel('$x$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
colorbar;
caxis([min_value max_value])

ax = gca;
ax.FontSize = FONT_SIZE;
set(gca,'fontname', FONT_NAME)

colormap jet
set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf, filename, 'Resolution',300)

filename = strcat('homo_', category, '_prediction_response.pdf');

fig2 = figure ('Position',[200 10 850 400]) ;  %680

subplot(2,3,1);
hold on;
title('0.32 s','interpreter','latex','FontSize',FONT_SIZE);
contourf(xm,ym,W4(6:21,6:36),50,'LineColor','none')
ylabel('$y$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
xlabel('$x$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
% colorbar;
caxis([min_value max_value])

ax = gca;
ax.FontSize = FONT_SIZE;
set(gca,'fontname', FONT_NAME)

subplot(2,3,2);
hold on;
title('0.33 s','interpreter','latex','FontSize',FONT_SIZE);
contourf(xm,ym,W5(6:21,6:36),50,'LineColor','none')
ylabel('$y$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
xlabel('$x$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
% colorbar;
caxis([min_value max_value])

ax = gca;
ax.FontSize = FONT_SIZE;
set(gca,'fontname', FONT_NAME)

subplot(2,3,3);
hold on;
title('0.34 s','interpreter','latex','FontSize',FONT_SIZE);
contourf(xm,ym,W6(6:21,6:36),50,'LineColor','none')
ylabel('$y$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
xlabel('$x$ [m]','Interpreter','Latex','FontSize',FONT_SIZE);
colorbar;
caxis([min_value max_value])

ax = gca;
ax.FontSize = FONT_SIZE;
set(gca,'fontname', FONT_NAME)

colormap jet

set(gcf, 'Position', get(0, 'Screensize'));
exportgraphics(gcf, filename, 'Resolution',300)