function [F_eff_ge,F_eff_b]  = force_eff_on_gamma (F_eff_history2)

global n_b; global n_ge; global size_of_b;

n_b2 = flip(n_b(1:(size_of_b(2)+1)));
n_b2 = [n_b2,n_b((size_of_b(2)+2):end)];

n_ge2 = flip(n_ge(1:(size_of_b(2)+2)));
n_ge2 = [n_ge2,n_ge((size_of_b(2)+3):end)];


F_eff_b = F_eff_history2(n_b2,:);

F_eff_ge = F_eff_history2(n_ge2,:);