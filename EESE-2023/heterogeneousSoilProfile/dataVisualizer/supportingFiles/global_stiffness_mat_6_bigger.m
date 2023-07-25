function [GK] = global_stiffness_mat_6_bigger (MM)

global wm; global hm; global h; global rho; global G;
global T; global freq; global ne; global nn; global nb;
global time_intg_type; global dt; global beta; global gamma; 
global s_loc; global Ntimestep;
global GM; global Keff; global GM_inv; 
global Le1; global Le2; global Le3; global Le4; global Le5; global Le6;


% Stiffness Matrix
Ke = [0.6667   -0.1667   -0.3333   -0.1667;
     -0.1667    0.6667   -0.1667   -0.3333;
     -0.3333   -0.1667    0.6667   -0.1667;
     -0.1667   -0.3333   -0.1667    0.6667];
Ke = Ke;%  * 4;% (h^2);
 
% Build Global Stiffness Matrix
GK = sparse(zeros(nn,nn));

% GK = zeros(nn,nn);
% Heterogeneous Domain
% Location 
L1 = [0,wm;0,hm];
L2 = [0,wm;30,60];
L3 = [0,wm;60,75];
L4 = [0,wm;75,90];
L5 = [0,wm;90,105];
L6 = [0,wm;105,120];

%%   Building the elements
Le1 = [L1(1,1)/h+1+(wm/h)*L1(2,1)/h:(L1(1,2)-h)/h+1+(wm/h)*L1(2,1)/h];
for i = 2:(L1(2,2)-L1(2,1))/h
  Le1(i,:) = Le1(i-1,:) + (wm/h);  
end
Le1 = Le1(:);
Le1 = sort(Le1);

Le2 = [L2(1,1)/h+1+(wm/h)*L2(2,1)/h:(L2(1,2)-h)/h+1+(wm/h)*L2(2,1)/h];
for i = 2:(L2(2,2)-L2(2,1))/h
  Le2(i,:) = Le2(i-1,:) + (wm/h);  
end
Le2 = Le2(:); Le2 = sort(Le2);

Le3 = [L3(1,1)/h+1+(wm/h)*L3(2,1)/h:(L3(1,2)-h)/h+1+(wm/h)*L3(2,1)/h];
for i = 2:(L3(2,2)-L3(2,1))/h
  Le3(i,:) = Le3(i-1,:) + (wm/h);  
end
Le3 = Le3(:); Le3 = sort(Le3);

Le4 = [L4(1,1)/h+1+(wm/h)*L4(2,1)/h:(L4(1,2)-h)/h+1+(wm/h)*L4(2,1)/h];
for i = 2:(L4(2,2)-L4(2,1))/h
  Le4(i,:) = Le4(i-1,:) + (wm/h);  
end
Le4 = Le4(:); Le4 = sort(Le4);

Le5 = [L5(1,1)/h+1+(wm/h)*L5(2,1)/h:(L5(1,2)-h)/h+1+(wm/h)*L5(2,1)/h];
for i = 2:(L5(2,2)-L5(2,1))/h
  Le5(i,:) = Le5(i-1,:) + (wm/h);  
end
Le5 = Le5(:); Le5 = sort(Le5);

Le6 = [L6(1,1)/h+1+(wm/h)*L6(2,1)/h:(L6(1,2)-h)/h+1+(wm/h)*L6(2,1)/h];
for i = 2:(L6(2,2)-L6(2,1))/h
  Le6(i,:) = Le6(i-1,:) + (wm/h);  
end
Le6 = Le6(:); Le6 = sort(Le6);

%   Taking out the same elements
for i = length(Le1):-1:1
    Lia = ismember(Le1(i),Le2);
    if Lia == 1
        Le1(i) = [];
    end
end

for i = length(Le1):-1:1
    Lia = ismember(Le1(i),Le3);
    if Lia == 1
        Le1(i) = [];
    end
end

for i = length(Le1):-1:1
    Lia = ismember(Le1(i),Le4);
    if Lia == 1
        Le1(i) = [];
    end
end

for i = length(Le1):-1:1
    Lia = ismember(Le1(i),Le5);
    if Lia == 1
        Le1(i) = [];
    end
end

for i = length(Le1):-1:1
    Lia = ismember(Le1(i),Le6);
    if Lia == 1
        Le1(i) = [];
    end
end
%% 
% Building the matrix
for k = 1:length(Le1)
    e = Le1(k);
    for i = 1:4
    for j = 1:4
        mi = MM(e,i); 
        mj = MM(e,j);
        GK(mi,mj) = Ke(i,j)*G(1)+ GK(mi,mj);
    end
    end
end 

for k = 1:length(Le2)
    e = Le2(k);
    for i = 1:4
    for j = 1:4
        mi = MM(e,i); 
        mj = MM(e,j);
        GK(mi,mj) = Ke(i,j)*G(2)+ GK(mi,mj);
    end
    end
end 

for k = 1:length(Le3)
    e = Le3(k);
    for i = 1:4
    for j = 1:4
        mi = MM(e,i); 
        mj = MM(e,j);
        GK(mi,mj) = Ke(i,j)*G(3)+ GK(mi,mj);
    end
    end
end 

for k = 1:length(Le4)
    e = Le4(k);
    for i = 1:4
    for j = 1:4
        mi = MM(e,i); 
        mj = MM(e,j);
        GK(mi,mj) = Ke(i,j)*G(4)+ GK(mi,mj);
    end
    end
end 

for k = 1:length(Le5)
    e = Le5(k);
    for i = 1:4
    for j = 1:4
        mi = MM(e,i); 
        mj = MM(e,j);
        GK(mi,mj) = Ke(i,j)*G(5)+ GK(mi,mj);
    end
    end
end 

for k = 1:length(Le6)
    e = Le6(k);
    for i = 1:4
    for j = 1:4
        mi = MM(e,i); 
        mj = MM(e,j);
        GK(mi,mj) = Ke(i,j)*G(6)+ GK(mi,mj);
    end
    end
end 