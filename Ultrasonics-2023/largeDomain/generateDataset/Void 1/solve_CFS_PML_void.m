function [u_full_history] = solve_CFS_PML_void (GK5,GC5,GM5,GF_m5,GF_m6)
%GK5 = GK; GC5 = GC; GM5 = GM; GF_m5 = GF_m;
global EBC_nn; global Ntimestep; global dt; global N_Node_aux; global NH; global NV;
global x_arr; global y_arr; global Lpml; global L;
global N_void; global xRotated_total; global yRotated_total;

%% Step 4 Tailor
total_nn = [1:length(GK5)];
All_nn = [1:length(GK5)];

for i = length(EBC_nn):-1:1
    All_nn(:,EBC_nn(1,i)) =[];
end

% Tailor GK, GC, GM, and GF
for i = length(EBC_nn):-1:1
    GK5(EBC_nn(1,i),:) =[];  % Delete row
    GK5(:,EBC_nn(1,i)) =[];  % Delete column
    
    GC5(EBC_nn(1,i),:) =[];  % Delete row
    GC5(:,EBC_nn(1,i)) =[];  % Delete column
    
    GM5(EBC_nn(1,i),:) =[];  % Delete row
    GM5(:,EBC_nn(1,i)) =[];  % Delete column
    
    GF_m5(EBC_nn(1,i),:) =[];
    GF_m6(EBC_nn(1,i),:) =[];
end


%% Step 5 : Solve u using bigger domain and GF_m
% [u_curt_history_big] = solve_u_DRM (GF_m);
load('MyColormap.mat')
ts = 0;
GF_t = GF_m5(:,ts+1);

% Keff = GM5 + 0.25 * GK5 * dt^2 + 0.5 * GC5 * dt;

u_curt = zeros(length(GF_t), 1);
udot_curt = zeros(length(GF_t), 1);

u_curt_history(:,1) = u_curt;

u_full = zeros(length(total_nn),1);
u_full(All_nn,1) = u_curt;
u_full_history(:,1) = u_full;

matrix_0= zeros(size(GK5));
matrix_I= eye(size(GK5));

J= -[matrix_0 -matrix_I; (GM5\GK5) matrix_0];
vector_0 = zeros(size(GF_t));
GF_k = (GF_m6(:, [1,2:2:end]));

%%
for ts = 1:Ntimestep
    
    u_prev = u_curt;
    udot_prev = udot_curt;
    
    s_prev = [u_prev;udot_prev];
    F1 = [vector_0;((GM5)\GF_m5(:,ts+1))];
    k1 = J*s_prev + F1;
    
    s_k2 = s_prev+0.5*dt*k1;
    F2 = [vector_0;((GM5)\GF_k(:,ts+1))];
    k2 = J*s_k2 + F2;
    
    s_k3 = s_prev+0.5*dt*k2;
    F3 = [vector_0;((GM5)\GF_k(:,ts+1))];
    k3 = J*s_k3 + F3;
    
    s_k4 = s_prev+dt*k3;
    F4 = [vector_0;((GM5)\GF_m5(:,ts+2))];
    k4 = J*s_k4 + F4;
    
    k = 1/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    s_curt = s_prev + dt*k;
    
    u_curt = s_curt(1:size(u_curt));
    udot_curt = s_curt(size(u_curt)+1:end);
    
    u_curt_history(:,ts+1) = u_curt;
    
    % Converting the solution vector to the solution matrix.
    u_full = sparse(zeros(length(total_nn),1));
    u_full(All_nn,1) = u_curt;
    u_full_history(:,ts+1) = u_full;
    
    
    if (mod(ts,10)==0)
        
        N_DOF_x = (length(u_full)/2);% - (N_Node_aux*2)*3)/2;  %length(u_curt)/2;
        N_DOF_y = N_DOF_x;
        for i = 1:(2*NV+1)
            soln_matrix_x(:,i) = u_full([(i-1)*(2*NH+1)+1:i*(2*NH+1)],1);
            soln_matrix_y(:,i) = u_full([N_DOF_x+(i-1)*(2*NH+1)+1:N_DOF_x+i*(2*NH+1)],1);
        end
        
        [row ,col] = size(soln_matrix_x);
        for i = 1:row
            for j = 1:col
                u_amplitude (i,j) =  sqrt( soln_matrix_x(i,j)^2 + soln_matrix_y(i,j)^2  );
            end
        end
        
        
        fig_handle = figure('Position',[200 200 600 400]);
        hold on;
        %
        contourf(x_arr,y_arr,u_amplitude',100,'LineColor','none');
        colorbar;
        caxis([0*10^-8 2*10^-8])
        
        xlabel('$x$ [m]', 'interpreter', 'latex')
        ylabel('$y$ [m]', 'interpreter', 'latex')
        
        colormap jet
        
        axis equal
        
        for p= 1:N_void
            plot(xRotated_total(:,p),yRotated_total(:,p));
            fill(xRotated_total(:,p),yRotated_total(:,p),'w')
        end
        
        filename =  ['graphs_u/u_at_' num2str(ts) '_th_iteration.pdf'];
        print(fig_handle, '-r90', '-dpdf', filename);
        hold off; clf; close all;
        
        
    end
    
end