function [GK,GM,GC] = building_matrices_hetero_void2 (Connectivity,Count_Inclusions)

global L; global H; global h; global NV; global NH;
global N_DOF; global N_Node; global N_Ele;
global rho; global G; global E; global nu;
global Vs; global Vp; global Mod_p; global Zs; global Zp;  global number_of_layers; global number_of_inclusions
global Li; global small_domain;
global ei; global Ele; global e_all; global thickness
global weight_total;
global EK_xx; global EK_xy; global EK_yx; global EK_yy; global EM_xx; global EM_yy;
global Layer_length;

%% Count inclusions
if Count_Inclusions == 'Y'
    inclusions = number_of_inclusions;
    for k = 1:inclusions
        L5 = cell2mat(Li(k));
        Le5 = [L5(1,1)/h+1+(L/h)*L5(2,1)/h:(L5(1,2)-h)/h+1+(L/h)*L5(2,1)/h];
        for i = 2:(L5(2,2)-L5(2,1))/h
            Le5(i,:) = Le5(i-1,:) + (L/h);
        end
        Le5 = Le5(:); Le5 = sort(Le5); ELi{k} = Le5';
    end
    
else
    inclusions = 0;
end

%% Plane-strain case
%% Orthotropic material
E_x1= 69*10^9;E_y1= 69*10^9; E_z1= 69*10^9; %% Material 01 (Aluminium)

v_xy1 = 0.33; v_yx1 = v_xy1*(E_y1/E_x1);
v_yz1 = 0.33; v_zy1 = v_yz1*(E_z1/E_y1);
v_xz1 = 0.33; v_zx1 = v_xz1*(E_z1/E_x1);
G_1 = (E_x1*E_y1)/(E_x1 + E_y1 +2*E_y1*v_xy1);

Q11 = (1-v_xz1*v_zx1)/E_x1;
Q12 = -(v_xy1+v_xz1*v_zy1)/E_x1;
Q13 = 0;
Q21 = -(v_yx1+v_zx1*v_yz1)/E_y1;
Q22 = (1- v_yz1*v_zy1)/E_y1;
Q23 = 0;
Q31 = 0;
Q32 = 0;
Q33 = 1/(G_1);

Q_matrix= [Q11 Q12 Q13; Q21 Q22 Q23; Q31 Q32 Q33];
C_matrix= inv(Q_matrix);

C11g(1) = C_matrix(1,1); % (E_x1*(v_xy1*v_yx1 - 1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C12g(1) = C_matrix(1,2); % -(E_y1*(v_xy1 + v_xz1*v_zy1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C21g(1) = C_matrix(2,1); % -(E_x1*(v_yx1 + v_yz1*v_zx1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C22g(1) = C_matrix(2,2); % (E_y1*(v_xz1^2 - 1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C33g(1) = C_matrix(3,3);

E_x2= 150*10^9;E_y2= 10*10^9; E_z2= 10*10^9; %% Material 02 (E-glass Epoxy)
v_xy2 = 0.33; v_yx2 = v_xy2*(E_y2/E_x2);
v_yz2 = 0.44; v_zy2 = v_yz2*(E_z2/E_y2);
v_xz2 = 0.33; v_zx2 = v_xz2*(E_z2/E_x2);
G_2 = (E_x2*E_y2)/(E_x2 + E_y2 +2*E_y2*v_xy2);

Q11 = (1-v_xz2*v_zx2)/E_x2;
Q12 = -(v_xy2+v_xz2*v_zy2)/E_x2;
Q13 = 0;
Q21 = -(v_yx2+v_zx2*v_yz2)/E_y2;
Q22 = (1- v_yz2*v_zy2)/E_y2;
Q23 = 0;
Q31 = 0;
Q32 = 0;
Q33 = 1/(G_2);

Q_matrix= [Q11 Q12 Q13; Q21 Q22 Q23; Q31 Q32 Q33];
C_matrix= inv(Q_matrix);

C11g(2) = C_matrix(1,1); % (E_x1*(v_xy1*v_yx1 - 1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C12g(2) = C_matrix(1,2); % -(E_y1*(v_xy1 + v_xz1*v_zy1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C21g(2) = C_matrix(2,1); % -(E_x1*(v_yx1 + v_yz1*v_zx1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C22g(2) = C_matrix(2,2); % (E_y1*(v_xz1^2 - 1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C33g(2) = C_matrix(3,3);

E_x3= 69*10^9;E_y3= 69*10^9; E_z3= 69*10^9; %% Material 03 (Aluminium)
v_xz3 = 0.33; v_zx3 = v_xz3*(E_z3/E_x3);
v_xy3 = 0.33; v_yx3 = v_xy3*(E_y3/E_x3);
v_yz3 = 0.33; v_zy3 = v_yz3*(E_z3/E_y3);
G_3 = (E_x3*E_y3)/(E_x3 + E_y3 +2*E_y3*v_xy3);

Q11 = (1-v_xz3*v_zx3)/E_x3;
Q12 = -(v_xy3*(v_yx3 + v_zx3*v_yz3))/(v_yx3* E_x3);
Q13 = 0;
Q21 = -(v_xy3*(v_yx3 + v_zx3*v_yz3))/(v_yx3* E_x3);
Q22 = (1- v_yz3*v_zy3)/E_y3;
Q23 = 0;
Q31 = 0;
Q32 = 0;
Q33 = 1/(G_3);

Q_matrix= [Q11 Q12 Q13; Q21 Q22 Q23; Q31 Q32 Q33];
C_matrix= inv(Q_matrix);

C11g(3) = C_matrix(1,1); % (E_x1*(v_xy1*v_yx1 - 1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C12g(3) = C_matrix(1,2); % -(E_y1*(v_xy1 + v_xz1*v_zy1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C21g(3) = C_matrix(2,1); % -(E_x1*(v_yx1 + v_yz1*v_zx1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C22g(3) = C_matrix(2,2); % (E_y1*(v_xz1^2 - 1))/(2*v_xy1*v_yx1 + v_xz1^2 - v_xy1*v_xz1^2*v_yx1 + v_xy1*v_yz1*v_zx1 + v_xz1*v_yx1*v_zy1 + v_xz1*v_yz1*v_zx1*v_zy1 - 1);
C33g(3) = C_matrix(3,3);

%% Step 2: Building Matrices

psi_xi(1) = {@(xi,eta) 1/4.*(2.*xi-1).*(eta.*eta-eta)};
psi_xi(2) = {@(xi,eta) 1/2.*(-2.*xi) .*(eta.*eta-eta)};
psi_xi(3) = {@(xi,eta) 1/4.*(2.*xi+1).*(eta.*eta-eta)};
psi_xi(4) = {@(xi,eta) 1/2.*(2.*xi-1).*(1-eta.*eta)};
psi_xi(5) = {@(xi,eta)      (-2.*xi) .*(1-eta.*eta)};
psi_xi(6) = {@(xi,eta) 1/2.*(2.*xi+1).*(1-eta.*eta)};
psi_xi(7) = {@(xi,eta) 1/4.*(2.*xi-1).*(eta.*eta+eta)};
psi_xi(8) = {@(xi,eta) 1/2.*(-2.*xi) .*(eta.*eta+eta)};
psi_xi(9) = {@(xi,eta) 1/4.*(2.*xi+1).*(eta.*eta+eta)};

psi_eta(1) = {@(xi,eta) 1/4.*(xi.*xi-xi).*(2.*eta-1)};
psi_eta(2) = {@(xi,eta) 1/2.*(1-xi.*xi) .*(2.*eta-1)};
psi_eta(3) = {@(xi,eta) 1/4.*(xi.*xi+xi).*(2.*eta-1)};
psi_eta(4) = {@(xi,eta) 1/2.*(xi.*xi-xi).*(-2.*eta)};
psi_eta(5) = {@(xi,eta)       (1-xi.*xi).*(-2.*eta)};
psi_eta(6) = {@(xi,eta) 1/2.*(xi.*xi+xi).*(-2.*eta)};
psi_eta(7) = {@(xi,eta) 1/4.*(xi.*xi-xi).*(2.*eta+1)};
psi_eta(8) = {@(xi,eta) 1/2.*(1-xi.*xi) .*(2.*eta+1)};
psi_eta(9) = {@(xi,eta) 1/4.*(xi.*xi+xi).*(2.*eta+1)};

psi(1) = {@(xi,eta) 1/4.*(xi.*xi-xi).*(eta.*eta-eta)};
psi(2) = {@(xi,eta) 1/2.*(1-xi.*xi) .*(eta.*eta-eta)};
psi(3) = {@(xi,eta) 1/4.*(xi.*xi+xi).*(eta.*eta-eta)};
psi(4) = {@(xi,eta) 1/2.*(xi.*xi-xi).*(1-eta.*eta)};
psi(5) = {@(xi,eta)      (1-xi.*xi) .*(1-eta.*eta)};
psi(6) = {@(xi,eta) 1/2.*(xi.*xi+xi).*(1-eta.*eta)};
psi(7) = {@(xi,eta) 1/4.*(xi.*xi-xi).*(eta.*eta+eta)};
psi(8) = {@(xi,eta) 1/2.*(1-xi.*xi) .*(eta.*eta+eta)}; %%%
psi(9) = {@(xi,eta) 1/4.*(xi.*xi+xi).*(eta.*eta+eta)};

%% Sub-dividing the layers
for k = 1:number_of_layers
    Ele{k} = (1+((NH*NV)*(Layer_length(k)/H))):((NH*NV)*(Layer_length(k+1)/H));
end

%% Step 3: Assembly of Mass Matrix

GK = sparse(N_DOF,N_DOF);
GM = sparse(N_DOF,N_DOF);

for layer = 1:(number_of_layers+inclusions) % Number of layers
    layer;
    EM_xx_1 = zeros(9,9);
    rho_layer= rho(layer);
    
    
    W= [1/3 4/3 1/3];
    xi= [-1 0 1];
    eta= [-1 0 1];
    for k= 1:3
        for l= 1:3
            for i = 1:9
                for j = 1:9
                    % For M matrix
                    value = W(k)*W(l)*(h/2)^2 * (psi{i}(xi(k),eta(l)).*psi{j}(xi(k),eta(l)));
                    EM_xx_1 (i,j) = EM_xx_1 (i,j) + rho_layer* value;
                    
                    EM_xx{1, layer} = EM_xx_1;
                    
                    
                end
            end
        end
    end
    EM_yy = EM_xx;
    
    Le = cell2mat(Ele(layer));
    for k = 1:length(Le)
        e = Le(k);
        
        for local_row = 1:9
            for local_col = 1:9
                % Assembly of GK
                % From EM_xx
                global_row = Connectivity(e,local_row );
                global_col = Connectivity(e,local_col );
                GM(global_row,global_col) = GM(global_row,global_col) + EM_xx{1, layer}(local_row,local_col);
                
                % From EM_xy
                global_row = Connectivity(e,local_row );
                global_col = Connectivity(e,local_col )+ N_Node ;
                GM(global_row,global_col) = GM(global_row,global_col) + 0;
                
                % From EM_yx
                global_row = Connectivity(e,local_row )+ N_Node;
                global_col = Connectivity(e,local_col );
                GM(global_row,global_col) = GM(global_row,global_col) + 0;
                
                % From EM_yy
                global_row = Connectivity(e,local_row)+ N_Node;
                global_col = Connectivity(e,local_col)+ N_Node;
                GM(global_row,global_col) = GM(global_row,global_col) + EM_yy{1, layer}(local_row,local_col);
                
            end
        end
    end
    
    
end

%%



if Count_Inclusions == 'Y'
    for jj = 1:inclusions
        for kk = 1:number_of_layers
            for i = length(Ele{kk}):-1:1
                Lia = ismember(Ele{kk}(i),ELi{jj});
                if Lia == 1
                    Ele{kk}(i) = [];
                end
            end
        end
    end
    Ele = [Ele,ELi];
end


%% Step 3 Assembly of K Matrix
for layer = 1:(number_of_layers+inclusions) %number_of_layers
    layer;
    C11 = C11g(layer);
    C12 = C12g(layer);
    C21 = C12;
    C22 = C22g(layer);
    C33 = C33g(layer);
    
    EK_xx1= zeros(9,9);
    EK_xx2= zeros(9,9);
    
    EK_yy1= zeros(9,9);
    EK_yy2= zeros(9,9);
    
    EK_xy1= zeros(9,9);
    EK_xy2= zeros(9,9);
    
    EK_yx1= zeros(9,9);
    EK_yx2= zeros(9,9);
    
    for k= 1:3
        for l= 1:3
            for i = 1:9
                for j = 1:9
                    value_1 = W(k)* W(l)*  psi_xi{i}(xi(k),eta(l)).*psi_xi{j}(xi(k),eta(l));
                    EK_xx1 (i,j) = EK_xx1(i,j) + C11*value_1;
                    EK_yy2 (i,j) = EK_yy2(i,j) + C33*value_1;
                    
                    value_2 = W(k) * W(l) * psi_eta{i}(xi(k),eta(l)).*psi_eta{j}(xi(k),eta(l));
                    EK_xx2 (i,j) = EK_xx2 (i,j) + C33*value_2;
                    EK_yy1 (i,j) = EK_yy1 (i,j) + C22*value_2;
                    
                    value_3 = W(k)* W(l)* psi_xi{i}(xi(k),eta(l)).*psi_eta{j}(xi(k),eta(l));
                    EK_xy1 (i,j) = EK_xy1 (i,j) + C12*value_3;
                    EK_yx2 (i,j) = EK_yx2 (i,j) + C33*value_3;
                    
                    value_4 = W(k)* W(l)* psi_eta{i}(xi(k),eta(l)).*psi_xi{j}(xi(k),eta(l));
                    EK_xy2 (i,j) = EK_xy2 (i,j) + C33*value_4;
                    EK_yx1 (i,j) = EK_yx1 (i,j) + C12*value_4;
                    
                    EK_xx{1, layer} = EK_xx1 + EK_xx2;
                    EK_xy{1, layer} = EK_xy1 + EK_xy2;
                    EK_yx{1, layer} = EK_yx1 + EK_yx2;
                    EK_yy{1, layer} = EK_yy1 + EK_yy2;
                    
                end
            end
        end
    end
    
    Le = cell2mat(Ele(layer));
    for k = 1:length(Le)
        e = Le(k);
        
        for local_row = 1:9
            for local_col = 1:9
                
                % From EK_xx
                global_row = Connectivity(e,local_row );
                global_col = Connectivity(e,local_col );
                GK(global_row,global_col) = GK(global_row,global_col) + EK_xx{1, layer}(local_row,local_col);
                
                % From EK_xy
                global_row = Connectivity(e,local_row );
                global_col = Connectivity(e,local_col )+ N_Node ;
                GK(global_row,global_col) = GK(global_row,global_col) + EK_xy{1, layer}(local_row,local_col);
                
                % From EK_yx
                global_row = Connectivity(e,local_row )+ N_Node;
                global_col = Connectivity(e,local_col );
                GK(global_row,global_col) = GK(global_row,global_col) + EK_yx{1, layer}(local_row,local_col);
                
                % From EK_yy
                global_row = Connectivity(e,local_row)+ N_Node;
                global_col = Connectivity(e,local_col)+ N_Node;
                GK(global_row,global_col) = GK(global_row,global_col) + EK_yy{1, layer}(local_row,local_col);
            end
            
        end
        
    end
    
end

%%
GC = sparse(N_DOF,N_DOF);