% 2D Plane-Strain and DRM  % Written by Bruno Guidio 
%Modified by Fazle Mahdi PRanto

clear all; clc;
global L; global H; global h; global NV; global NH; global rho; global G; global E; global nu; 
global Vs; global Vp; global Mod_p; global Zs; global Zp;  
global N_Ele; global N_Node; global N_DOF; global dt; global Ntimestep; 
global GC; global x_arr; global y_arr; global EBC_nn; global s_loc; global um; 
global number_of_layers; global number_of_inclusions
global e_all; global e_pml; global ei; global Lpml; global Ele; global Li; global N_Node_aux;
 global xy_coord_node; global N_iteration;
global N_void; global Eli_fea; global xRotated_total; global yRotated_total;
global weight_total; global Crack_data_yes_no_1; 
global EK_xx; global EK_xy; global EK_yx; global EK_yy; global EM_xx; global EM_yy;
global Layer_length;
mkdir('graphs_u'); 
% mkdir('graphs_ux')
% mkdir('graphs_uy')
mkdir('graphs_void')
load('MyColormap.mat')
LASTN = maxNumCompThreads(1);


% Step 1: Dimensions

%size_of_interior_domain = [0.6,0.7];                  %[x,y]size of domain without PML
%Lpml = 0.1;

L = 0.1;               %Total Length of the domain
H = 0.05;               %Total Height of the domain

h = 0.001;                %Element size 
Layer_length= [0, 0.015, 0.035, 0.05];  %Start with zero and then cummulative.
 

number_of_layers = 3;                               %In my current code, my layers can't have different sizes, 
rng('shuffle')                                      %Random number generator setting                                                    %so we have to use: number_of_layers = (H/Lpml)

number_of_inclusions = 0;                           %if we want to have inclusions. Currently, this is 0.                                   

dt = 5*10^-8;                                       
Ntimestep = 350;                           %Number of Timesteps

% >>>>>>>> <<<<<<<<<
Ricker_freq = 500000;                                %Hz {frequency - we can go up to 20 Hz}
point_load_loc_1 = [0, 0.025];                         %Location of the force [x,y] - 





                                                    %the x and y coordinates are on the bottom-left corner of the domain

sensor_x = [0.01:0.002:0.09]';
sensor_y(1:41, 1) = 0.05;

sensor_loc= [sensor_x, sensor_y];



NV = single(H/h);                                   %Number of elements in the vertical axis
NH = L/h;                                           %Number of elements in the horizontal axis 
N_Ele = NV * NH;                                    %Total number of elements
N_Node = (2*NV+1)*(2*NH+1);                         %Total Number of nodes
N_DOF = 2*N_Node;                                   %Number of degrees of freedom before PML

%Loading information 


Point_loading_x_1 = point_loading_1(Ricker_freq);       %Function to create the loading 
Point_loading_x_2 = point_loading_2(Ricker_freq);
    

% For postprocessing. 
x_arr = [0:1:2*(NH)];
x_arr = x_arr *h/2;

y_arr = [0:1:2*(NV)];
y_arr = y_arr *h/2;

%% Material property

%Soil
rho = [2700, 1550, 2700];                                       %Mass density of the entire domain
nu = 0.33;
E = 69*10^9;   
 
G = E./(2.*(1+nu));
Mod_p = 2*G.*(1 - nu)./(1 - 2*nu);  
Vs = sqrt(G./rho);
Vp = sqrt(Mod_p./rho);
Zp = rho.*Vp;
Zs = rho.*Vs;   

%% Enter connectivity array.
[Connectivity,EBC_nn] = connectivity_matrix ();     %Connectivity array linking every element with its corresponding nodes
%% Defining what elements are in the PML layer

e_all = 1:(N_Ele);


%% Matrices
Count_Inclusions = 'N';
[GK,GM,GC] = building_matrices_hetero_void2 (Connectivity,Count_Inclusions);


%% RANDOMIZER
N_iteration = 2000;                                    %Iterations in the randomizer1
Crack_data_yes_no_1= zeros(N_iteration,N_Ele);      %Save location of crack in each iteration

for iter = 1:N_iteration                            %For loop for randomizer
%% Void features:
tic
save_void = [];
GM_new= [];
GK_new= [];

void_features_test();                               %function to define void features
[weight_total_per_element] = void_elements (Connectivity,iter);                  %function to define each elements are void-elements
count = 1;

if N_void >= 1
    for i= 1: (NV*NH)
        if weight_total_per_element(1,i)<= 0.0001
            save_void(count,1)= i;
            count= count +1;
        end
    end
else
    save_void(1,1)= 0;
end


%%Converting void element layer by layer
% for element = 1: length(save_void)
%     if save_void(element,1)<= length(Ele{1,1})
%         temp_lay_1(element,1)= save_void(element,1);
%         
%     elseif save_void(element,1)> length(Ele{1,1}) && save_void(element,1)<= length(Ele{1,2})
%         temp_lay_1(element,2)= save_void(element,1);
%     else 
%         temp_lay_1(element,3)= save_void(element,1);
%     end
%     
% end




%%

GM_new = GM ;
GK_new = GK ;

%% building Mass void matrices:

if save_void== 0
    GM_new = GM ;
else

    
    for k = 1: length(save_void)  %%Figure out how to automate this.
        e = save_void(k,1);
        if e<= length(Ele{1,1})
        EM_xx_layer = EM_xx {1,1};
        EM_yy_layer = EM_yy {1,1};

        elseif e> length(Ele{1,1}) && e<= length(Ele{1,1})+ length(Ele{1,2})
        EM_xx_layer = EM_xx {1,2};
        EM_yy_layer = EM_yy {1,2};

        else 
        EM_xx_layer = EM_xx {1,3};
        EM_yy_layer = EM_yy {1,3};

        end
        for local_row = 1:9
            for local_col = 1:9
                global_row = Connectivity(e,local_row );
                global_col = Connectivity(e,local_col );
                GM_new(global_row, global_col)= (GM_new(global_row, global_col)-EM_xx_layer(local_row, local_col))+ EM_xx_layer(local_row, local_col)* weight_total_per_element(e);
            
            
                global_row = Connectivity(e,local_row )+ N_Node;
                global_col = Connectivity(e,local_col )+ N_Node;
                GM_new(global_row, global_col)=(GM_new(global_row, global_col)-EM_yy_layer(local_row, local_col))+ EM_yy_layer(local_row, local_col)* weight_total_per_element(e);
           
            end
        end
    end

end

%% building Stiffness void matrices:
if save_void== 0
    GK_new = GK;
else
for k = 1: length(save_void)
    e = save_void(k);
    if e<= length(Ele{1,1})
        EK_xx_layer = EK_xx {1,1};
        EK_xy_layer = EK_xy {1,1};
        EK_yx_layer = EK_yx {1,1};
        EK_yy_layer = EK_yy {1,1};
    elseif e> length(Ele{1,1}) && e<= length(Ele{1,1}) + length(Ele{1,2})
        EK_xx_layer = EK_xx {1,2};
        EK_xy_layer = EK_xy {1,2};
        EK_yx_layer = EK_yx {1,2};
        EK_yy_layer = EK_yy {1,2};
    else 
        EK_xx_layer = EK_xx {1,3};
        EK_xy_layer = EK_xy {1,3};
        EK_yx_layer = EK_yx {1,3};
        EK_yy_layer = EK_yy {1,3};
    end
    for local_row = 1:9
        for local_col = 1:9
            global_row = Connectivity(e,local_row );
            global_col = Connectivity(e,local_col );
            GK_new(global_row, global_col)= (GK_new(global_row, global_col)-EK_xx_layer(local_row, local_col))+ EK_xx_layer(local_row, local_col)* weight_total_per_element(e);
            
            
            global_row = Connectivity(e,local_row );
            global_col = Connectivity(e,local_col )+ N_Node;
            GK_new(global_row, global_col)= (GK_new(global_row, global_col)-EK_xy_layer(local_row, local_col))+ EK_xy_layer(local_row, local_col)*weight_total_per_element(e);
            
            
            global_row = Connectivity(e,local_row )+ N_Node;
            global_col = Connectivity(e,local_col );
            GK_new(global_row, global_col)= (GK_new(global_row, global_col)-EK_yx_layer(local_row, local_col))+ EK_yx_layer(local_row, local_col)* weight_total_per_element(e);
            
            
            global_row = Connectivity(e,local_row )+ N_Node;
            global_col = Connectivity(e,local_col )+ N_Node;
            GK_new(global_row, global_col)= (GK_new(global_row, global_col)-EK_yy_layer(local_row, local_col))+  EK_yy_layer(local_row, local_col)* weight_total_per_element(e);
            
        end
    end
end
end
        

%% Body force
GF4 = sparse(zeros(length(GK),1));
GF_m = zeros(length(GK),Ntimestep+2);

GF_m((2*NH+1)*(2*point_load_loc_1(2)/h)+1+(2*point_load_loc_1(1)/h),:) = ...
                      GF4((2*NH+1)*(2*point_load_loc_1(2)/h)+1+(2*point_load_loc_1(1)/h),:) - Point_loading_x_1';
                  

                  
GF_m2 = zeros(length(GK),2*Ntimestep+1);

GF_m2((2*NH+1)*(2*point_load_loc_1(2)/h)+1+(2*point_load_loc_1(1)/h),:) = ...
                      GF4((2*NH+1)*(2*point_load_loc_1(2)/h)+1+(2*point_load_loc_1(1)/h),:) - Point_loading_x_2';
                  

%% Solving the forward problem
[u_full_history_big] = solve_CFS_PML_void (GK_new,GC,GM_new,GF_m,GF_m2);

%% Saving ux for each sensor
um = zeros(length(sensor_loc),Ntimestep+1);       %um
for i = 1:length(sensor_loc)
    um(i,:) = u_full_history_big(round(N_Node+(2*NH+1)*(2*sensor_loc(i,2)/h)+1+(2*sensor_loc(i,1)/h)),:);
% um(i+length(sensor_loc),:) = u_full_history_big(round(N_Node+(2*NH+1)*(2*sensor_loc(i,2)/h)+1+(2*sensor_loc(i,1)/h)),:);
end
%%

%% Code to chek if [um] has any NAN -- if YES, the code will stop
% Pranto can comment these lines
ans_nan = isnan(um);
ans_nan = double(ans_nan);

k = find(ans_nan);

if length(k) > 0
    break
end
%%
um_history{iter,1}= um;

Eli_fea_history{iter,1} = Eli_fea;              %if needed


test_avg= nonzeros(Crack_data_yes_no_1(iter,:));     
weight_avg= length(test_avg);                        
Crack_per_iteration(iter,1)= weight_avg;             

% fig_handle2 = figure('Position',[200 200 600 400]);
% hold on; grid on;
% %Creating Void:
%     for p= 1:N_void
%         plot(xRotated_total(:,p),yRotated_total(:,p));
%         fill(xRotated_total(:,p),yRotated_total(:,p),'r')
%     end
%         
% xlabel('x[m]')
% ylabel('y[m]')
% xlim([0 L])
% ylim([0 H])  
%         
% filename =  ['graphs_void/iteration_num_' num2str(iter) '.pdf']; 
% print(fig_handle2, '-r90', '-dpdf', filename);
% hold off; clf; close all;

toc
iter
end                                     %end of loop iteration 

%% After all the iterations - let's save the data
Average_crack_element = round(mean(Crack_per_iteration));
for x= 1: N_Ele
    test= nonzeros(Crack_data_yes_no_1(:,x));
    weight_total= length(test);
    Crack_data_total(1,x)= weight_total;
end
crack_percentage= round((Average_crack_element/N_Ele)*100)
 save('Crack_data_void_0.mat','Crack_data_yes_no_1')
 save('Signal_ani_void_0.mat','um_history')
%  save('Eli_fea_all.mat','Eli_fea_history')

 %%
  cdata=[];
for i=1:(NV)
    for j=1:(NH)
        e = (i-1)*NH +j;
        cdata(i,j)= Crack_data_total (e);
    end
end

cdata= flipud(cdata);
mainDataAxes = axes('position',[0.08  0.08  0.85    0.85]);
axes(mainDataAxes);
h_3 = imagesc(cdata);
J= jet;
J(1,:)= [1 1 1];
colormap(J);
colorbar;
set(gca, 'XTick',[20],'XTickLabel',{'X axis\rightarrow'},...
    'TickLength', [0 0]);
set(gca, 'YTick',[]);
ylabel('Y axis\rightarrow')
title(['Heat map'],'FontSize', 20);


set(gcf,'position',[100 100 800 500]);  % need to un-dock the figure window first
figure_size = get(gcf, 'position');
set(gcf,'PaperPosition',figure_size/100); 
 print(gcf,'-dpng','-r100', ['graphs_void/2D_heatmap_void_1.png']);
