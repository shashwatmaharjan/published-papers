function [F_eff_history2] = solve_F_eff_bigger (GF_m)

global wm; global hm; global h; global nn;
global corner_of_b; global size_of_b; global corner_of_i; global size_of_i;
global corner_of_e; global size_of_e;
global time_intg_type; global dt; global beta; global gamma; global Ntimestep;
global GK; global GM; global Meff; global GM_inv; global GC;
global n_i; global n_b; global n_e; global n_ge; global n_e_out; global n_matrix;

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

%%
n_small = n_matrix(corner_of_e(2)/h+1:(corner_of_e(2)/h+size_of_e(2)/h+1),...
    corner_of_e(1)/h+1:(corner_of_e(1)/h+size_of_e(1)/h+1));

n_small = reshape(n_small,1,(size_of_e(1)/h+1)*(size_of_e(2)/h+1));

%%  
Mbe = GM(n_b,n_e);
Kbe = GK(n_b,n_e);

Meb = GM(n_e,n_b);
Keb = GK(n_e,n_b);

%%
ts = 0;
GF = GF_m(:,ts+1);

% Newmark Time Integration
Meff = GM + 0.5 * GC * dt + 0.25 * GK * dt^2;

u_curt = sparse(zeros(length(GF),1));
udot_curt = sparse(zeros(length(GF),1));
uddot_curt = GM\GF;

Fb_eff = -Mbe*uddot_curt(n_e) -Kbe*u_curt(n_e);
Fe_eff =  Meb*uddot_curt(n_b) +Keb*u_curt(n_b);

F_eff = [zeros(length(n_i),1);Fb_eff;Fe_eff];

F_eff_history(:,1) = F_eff;

%%
for ts = 1: Ntimestep

    t = ts*dt;
    u_prev = u_curt;
    udot_prev = udot_curt;
    uddot_prev = uddot_curt;

    % Build the RHS load vector
    % Global Force Vector

    GF=  GF_m(:,ts+1);

    if (time_intg_type ==1)
        % Implicit time integration
        RHS = GF - GC * (udot_prev + uddot_prev * 0.5 * dt)  ...
            - GK*(u_prev + udot_prev*dt + 0.25*uddot_prev*dt^2);
        % Solve for acceleration at the current (i)-th time step
        uddot_curt = Meff\RHS;
    else
        % explicit time integration.
        RHS = GF - GK * (u_prev);
        uddot_curt = GM_inv * RHS;
    end

    % Updating u and udot at the current (i)-th time step
    u_curt = u_prev + udot_prev * dt + ...
        0.5 * (uddot_prev * (1.0-2.0*beta) + uddot_curt* 2.0*beta)*dt^2;

    udot_curt = udot_prev + ...
        (uddot_prev * (1.0-gamma) + uddot_curt* gamma)*dt;


    Fb_eff = -Mbe*uddot_curt(n_e) -Kbe*u_curt(n_e);
    Fe_eff =  Meb*uddot_curt(n_b) +Keb*u_curt(n_b);

    F_eff = [zeros(length(n_i),1);Fb_eff;Fe_eff];


    F_eff_history(:,ts+1) = F_eff;


end

%%
order_eff = [n_i';n_b';n_e'];

for i = 1:length(order_eff)
    k = find(order_eff==i);
    true_pos(i,1) = k;
end
true_pos = sparse(true_pos);

F_eff_history = F_eff_history(true_pos,:);

%%
F_eff_history2 = F_eff_history(n_small,:);

n_matrix2 = sparse([1:(size_of_e(1)/h+1)]);
for i = 2:(size_of_e(2)/h+1)
    n_matrix2(i,:) = n_matrix2(i-1,:)+n_matrix2(1,end);
end
n_small2 = reshape(n_matrix2,1,(size_of_e(1)/h+1)*(size_of_e(2)/h+1));

order_eff_small = [n_small2'];

for i = 1:length(order_eff_small)
    k = find(order_eff_small==i);
    true_pos2(i,1) = k;
end
true_pos2 = sparse(true_pos2);

F_eff_history2 = F_eff_history2(true_pos2,:);