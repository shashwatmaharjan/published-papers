function [GK, GM, GC] = FEM_matrices_bigger ( )

global wm; global hm; global h; global rho; global G;
global T; global freq; global ne; global nn; global nb;
global time_intg_type; global dt; global beta; global gamma;
global s_loc; global Ntimestep; global w_speed;
global Keff; global GM_inv; global MM; global Soil_profile; global Applied_force;

%% Mapping Matrix
nM = (hm/h)+1;  %number of mapping matrix
nRM = (wm/h);   %number of rows for each mapping matrix

MM = sparse(zeros((nM-1)*nRM,4)); %Sparse
MM(1,1:4) = [nRM+2, nRM+3, 2, 1];
for i = 2:nRM
    MM(i,1:4) = MM(i-1,1:4)+1;
end

for i = 2:nM-1
    MM(nRM*(i-1)+1:nRM*i,1:4) = MM(nRM*(i-2)+1:nRM*(i-1),1:4)+(wm/h+1);
end

GM = sparse((zeros(nn,nn))); %SPARSE

%%%%%%%%%%  Global K  %%%%%%%%%%
if (strcmp(Soil_profile,'Homogeneous'))
    GK = global_stiffness_homo(MM);
elseif (strcmp(Soil_profile,'Material_profile_6_bigger'))
    GK = global_stiffness_mat_6_bigger(MM);
elseif (strcmp(Soil_profile,'Material_profile_6'))
    GK = global_stiffness_mat_6(MM);
elseif (strcmp(Soil_profile,'Material_profile_8b_bigger'))
    GK = global_stiffness_mat_8b_bigger(MM);
elseif (strcmp(Soil_profile,'Material_profile_4'))
    GK = global_stiffness_mat_4(MM);
    % elseif (strcmp(Soil_profile,'Material_profile_5_2'))
    %     GK = global_stiffness_mat_5_2(MM);
    % elseif (strcmp(Soil_profile,'Material_profile_5_3'))
    %     GK = global_stiffness_mat_5_3(MM);
end

%% Mass Matrix
Me =  [0.4444    0.2222    0.1111    0.2222;
    0.2222    0.4444    0.2222    0.1111;
    0.1111    0.2222    0.4444    0.2222;
    0.2222    0.1111    0.2222    0.4444];
[row, col] = size(Me);

if (time_intg_type ==1)
    Me = rho * ((h^2)/4) * Me;
    Me = sparse(Me);
else %For Explicit
    Me_lumping = sparse(zeros(row,col));
    for i = 1:row
        for j = 1:col
            Me_lumping(i,i) = Me_lumping(i,i) + Me(i,j);
        end
    end
    Me = rho * ((h^2)/4) * Me_lumping;
end

%% Build Global Mass Matrix
for e = 1:ne
    for i = 1:4
        for j = 1:4
            mi = MM(e,i);
            mj = MM(e,j);
            GM(mi,mj) = Me(i,j)+ GM(mi,mj);
        end
    end
end

%% Global Mass Matrix - Inverse
[row, col] = size(GM);

if (time_intg_type ==1)
    GM_inv = 0;
else
    for i = 1:row
        GM_inv(i,i)= 1/GM(i,i);
    end
end

%% Global Dumping Matrix
if (strcmp(Soil_profile,'Material_profile_6_bigger'))
    GC = global_damping_abs_3_bigger(MM);
elseif (strcmp(Soil_profile,'Material_profile_8b_bigger'))
    GC = global_damping_abs_3_bigger(MM);
elseif (strcmp(Soil_profile,'Material_profile_6')) || (strcmp(Soil_profile,'Material_profile_4'))
    GC = global_damping_abs_3(MM);
    % elseif (strcmp(Soil_profile,'Material_profile_5')) ||  (strcmp(Soil_profile,'Material_profile_5_2'))||  (strcmp(Soil_profile,'Material_profile_5_3'))
    %     GC = global_damping_abs_3(MM);
else
    GC = global_damping_abs(MM);
    GC = (G(1)/w_speed(1)) * (h/2) * GC;
end
