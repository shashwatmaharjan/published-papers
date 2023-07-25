function [F_eff_ge,F_eff_b]  = force_eff_on_gamma (F_eff_history2)

global n_b; global n_ge; global size_of_b;

n_b2 = flip(n_b(1:(size_of_b(2)+1))); 
n_b2 = [n_b2,n_b((size_of_b(2)+2):end)];

n_ge2 = flip(n_ge(1:(size_of_b(2)+2)));
n_ge2 = [n_ge2,n_ge((size_of_b(2)+3):end)];


F_eff_b = F_eff_history2(n_b2,:);

F_eff_ge = F_eff_history2(n_ge2,:);


% %Grid Matrix
%     xb = [0:290];
%     xge = [0:294];
%     yb = [1:Ntimestep+1];
%     yge = [1:Ntimestep+1];
%     [xb,yb] = meshgrid(xb,yb);
%     [xge,yge] = meshgrid(xge,yge);
%    figure_F = figure(1); 
%    subplot(1,2,1)
%    hold on
%    %temp = [num2str(iteration_so_far),' -th iteration - TARGET']; 
%    %title(temp); 
%    contourf(xge,yge,F_eff_ge',50,'LineColor','none')
%    colorbar;
%    caxis([-1500 1500])
%    
%    subplot(1,2,2)
%    %F_g_plot = force_plot(F_g);
%    hold on
% %    temp = [num2str(iteration_so_far),' -th iteration - GUESS']; 
% %    title(temp); 
%    contourf(xb,yb,F_eff_b',50,'LineColor','none')
%    colorbar;
%    caxis([-1500 1500])