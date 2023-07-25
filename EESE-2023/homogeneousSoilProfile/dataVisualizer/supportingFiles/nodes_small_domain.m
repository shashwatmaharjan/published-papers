function [] = nodes_small_domain()

global wm; global hm; global h; global nn;
global corner_of_b; global size_of_b; global corner_of_i; global size_of_i;
global corner_of_e; global size_of_e;
global time_intg_type; global dt; global beta; global gamma; global Ntimestep;
global GK; global GM; global Meff; global GM_inv; global GC;
global n_i; global n_b; global n_e; global n_ge; global n_e_out;

%%
nn_i = (size_of_i(1)/h+1)*(size_of_i(2)/h+1);
nn_b = (size_of_b(1)/h+1)*(size_of_b(2)/h+1) - nn_i;
nn_e = nn - nn_b - nn_i;

n_matrix = sparse([1:(wm/h+1)]);
for i = 2:(hm/h+1)
    n_matrix(i,:) = n_matrix(i-1,:)+n_matrix(1,end);
end

n_i = n_matrix(corner_of_i(2)/h+1:(corner_of_i(2)/h+size_of_i(2)/h+1),...
    corner_of_i(1)/h+1:(corner_of_i(1)/h+size_of_i(1)/h+1));

n_i = reshape(n_i,1,nn_i);


n_b = n_matrix(corner_of_b(2)/h+1:(corner_of_b(2)/h+size_of_b(2)/h+1),...
    corner_of_b(1)/h+1:(corner_of_b(1)/h+size_of_b(1)/h+1));

n_b = reshape(n_b,1,nn_b+nn_i);

for i = length(n_b):-1:1
    Lia = ismember(n_b(i),n_i);
    if Lia == 1
        n_b(i) = [];
    end
end

n_e = reshape(n_matrix,1,nn);
for i = length(n_e):-1:1
    Lia = ismember(n_e(i),n_i);
    if Lia == 1
        n_e(i) = [];
    end
end
for i = length(n_e):-1:1
    Lia = ismember(n_e(i),n_b);
    if Lia == 1
        n_e(i) = [];
    end
end

%% Nodes on Gamma E (We will use these nodes when we solve the adj problem)
nn_ge = ((size_of_b(1)+2)/h+1)*((size_of_b(2)+1)/h+1) - nn_i - nn_b;

n_ge = n_matrix((corner_of_b(2)-1)/h+1:((corner_of_b(2)-1)/h+(size_of_b(2)+1)/h+1),...
    (corner_of_b(1)-1)/h+1:((corner_of_b(1)-1)/h+(size_of_b(1)+2)/h+1));

n_ge = reshape(n_ge,1,nn_ge+nn_b+nn_i);
for i = length(n_ge):-1:1
    Lia = ismember(n_ge(i),n_i);
    if Lia == 1
        n_ge(i) = [];
    end
end
for i = length(n_ge):-1:1
    Lia = ismember(n_ge(i),n_b);
    if Lia == 1
        n_ge(i) = [];
    end
end

n_e_out = n_e;
for i = length(n_e_out):-1:1
    Lia = ismember(n_e_out(i),n_ge);
    if Lia == 1
        n_e_out(i) = [];
    end
end