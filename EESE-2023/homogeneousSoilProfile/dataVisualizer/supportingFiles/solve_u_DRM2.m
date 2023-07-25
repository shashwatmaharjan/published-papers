function [u_curt_history] = solve_u_DRM2 (F_eff_history)

global wm; global hm; global h; global nn;
global corner_of_b; global size_of_b; global corner_of_i; global size_of_i;
global time_intg_type; global dt; global beta; global gamma; global Ntimestep;
global GK; global GM; global Meff; global GM_inv; global GC;
global n_i; global n_b; global n_e; global s_loc; global Meff_inv;

%%
ts = 0;
GF = F_eff_history(:,ts+1);

% Newmark Time Integration
Meff = GM + 0.5 * GC * dt + 0.25 * GK * dt^2;

u_curt = (zeros(length(GF),1));
udot_curt = (zeros(length(GF),1));
uddot_curt = GM\GF;

u_curt_history(:,1) = u_curt;

%%
for ts = 1: Ntimestep

    t = ts*dt;
    u_prev = u_curt;
    udot_prev = udot_curt;
    uddot_prev = uddot_curt;

    % Build the RHS load vector
    % Global Force Vector

    GF = F_eff_history(:,ts+1);
    % GF = GF(true_pos);

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

    u_curt_history(:,ts+1) = u_curt;

end