function [weight_total] = void_elements (Connectivity,iter)

global N_Ele; global xy_coord_node; global N_void; global Eli_fea;
global Crack_data_yes_no_1; global weight_total

for i = 1:N_Ele
    X_square(i,:)= xy_coord_node(Connectivity(i,:),1); % x-coordinate of each node of  element
    Y_square(i,:)= xy_coord_node(Connectivity(i,:),2); % y-coordinate of each node of  element
    
    if N_void == 0
        weight_total(i) = 1;
    else
        
        for p= 1:N_void
            
            for nnodes = 1:9
                result_square(nnodes)= isellipse(Eli_fea(p,3),Eli_fea(p,4),X_square(i,nnodes),Y_square(i,nnodes),Eli_fea(p,1),Eli_fea(p,2),Eli_fea(p,5));
                
                logic (nnodes) = result_square(nnodes)<1;
                
            end
            
            logic_center = logic(5);
            logic_others = logic([1:4,6:9]);
            others_inside = 0;
            for any_node = 1:8
                if logic_others(any_node) == 1
                    others_inside = others_inside + 1;
                end
            end
            
            if logic_center==1 && others_inside >= 4
                Crack_data_yes_no_1(iter,i)= 1;
                weight(1,p)= 0.0001;
            else
                weight(1,p)= 1;
            end
            
        end
        
        weight_total(i) = prod(weight);
        
    end
    
end