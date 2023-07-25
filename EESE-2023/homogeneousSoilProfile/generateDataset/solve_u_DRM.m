function [u_sensors] = solve_u_DRM(F_eff_history)

global wm; global hm; global h; global nn;
global corner_of_b; global size_of_b; global corner_of_i; global size_of_i;
global time_intg_type; global dt; global beta; global gamma; global Ntimestep;
global GK; global GM; global Meff; global GM_inv; global GC; global Meff_inv;
global n_i; global n_b; global n_e; global s_loc;

%%
ts = 0;
u_curt = (zeros(length(GM),1));
udot_curt = (zeros(length(GM),1));
uddot_curt = GM\F_eff_history(:,ts+1);

u_sensors (:,1) = u_curt(s_loc,1);

%%
for ts = 1: Ntimestep

    t = ts*dt;
    u_prev = u_curt;
    udot_prev = udot_curt;
    uddot_prev = uddot_curt;

    % Build the RHS load vector
    % Global Force Vector

    GF = F_eff_history(:,ts+1);
    RHS = GF - GC * (udot_prev + uddot_prev * 0.5 * dt)  ...
        - GK*(u_prev + udot_prev*dt + 0.25*uddot_prev*dt^2);

    uddot_curt = Meff\RHS;

    % Updating u and udot at the current (i)-th time step
    u_curt = u_prev + udot_prev * dt +  0.5 * (uddot_prev * (1.0-2.0*beta) + uddot_curt* 2.0*beta)*dt^2;

    udot_curt = udot_prev + (uddot_prev * (1.0-gamma) + uddot_curt* gamma)*dt;


    u_sensors (:,ts+1) = u_curt(s_loc,1);

end