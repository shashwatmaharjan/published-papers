function [Connectivity,EBC_nn] = connectivity_matrix ()
global NV; global NH; global N_Node; global EBC_nn;
global h; global xy_coord_ele; global xy_coord_node;

%% Step 1: Enter connectivity array.
EBC_nn1 =[]; EBC_nn2 =[]; Connectivity = [];
for i = 1:(NV)
    for j = 1:(NH)
        e = (i-1)*NH +j;
        x_coord_ele(e) = (j-1)*h + h/2 ;   %x-coordinate of center of each element
        y_coord_ele(e) = (i-1)*h + h/2 ;   %y-coordinate of center of each element
        if (j == 1)
            Connectivity(e,1) = (2*(i-1)) * (2*NH+1)+1;
            Connectivity(e,2) = Connectivity(e,1)+1; 
            Connectivity(e,3) = Connectivity(e,2)+1;           
            Connectivity(e,4) = (2*i-1)  * (2*NH+1)+1;
            Connectivity(e,5) = Connectivity(e,4)+1;
            Connectivity(e,6) = Connectivity(e,5)+1;            
            Connectivity(e,7) = (2*i) * (2*NH+1)+1;
            Connectivity(e,8) = Connectivity(e,7)+1;
            Connectivity(e,9) = Connectivity(e,8)+1;
            %Essential boundary condition array: when the LHS(left hand side) is fixed.
            if (i == 1) %the bottom row.
                EBC_nn1 = [EBC_nn1,[Connectivity(e,1),Connectivity(e,4),Connectivity(e,7)]];
            else
                EBC_nn1 = [EBC_nn1,[Connectivity(e,4),Connectivity(e,7)] ];
            end
        elseif (j==NH)
            for k = 1:9
                Connectivity(e,k) = Connectivity(e-1,k)+2;
            end
            if (i == 1) %the bottom row. %Right side fixed.
                EBC_nn2 = [EBC_nn2,[Connectivity(e,3),Connectivity(e,6),Connectivity(e,9)]];
            else
                EBC_nn2 = [EBC_nn2,[Connectivity(e,6),Connectivity(e,9)] ];
            end
        else
            for k = 1:9
                Connectivity(e,k) = Connectivity(e-1,k)+2;
            end
        end
    end
end

xy_coord_ele = [x_coord_ele; y_coord_ele]';

EBC_nn3 = [1:2*NH+1]; %Bottom side fixed.

% % EBC_nn1 = left side; EBC_nn2 = right side; EBC_nn3 = bottom
%EBC_nn = unique([EBC_nn1,EBC_nn2,EBC_nn3]);     %left, right, and bottom fixed!

EBC_nn = unique([EBC_nn2]);                     %only right side fixed


%% Double up the size of EBC_nn array (considering both x and y components of zero displacement)
temp = length(EBC_nn);
for i = 1:temp
    EBC_nn(1,i+temp) = EBC_nn(1,i)+N_Node;
end

%% Coordinates of each node
for i=1:2*(NV)+1
    for j=1:2*(NH)+1
        l= (i-1)*(2*NH+1) +j;
        x_coord_node(l) = (j-1)/2 *h;  %double-check if element size(h) is different of 1
        y_coord_node(l) = (i-1)/2 *h;
    end
end
xy_coord_node= [x_coord_node; y_coord_node]';