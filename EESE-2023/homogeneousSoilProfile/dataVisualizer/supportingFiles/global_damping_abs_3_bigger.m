function [GC] = global_damping_abs_3_bigger(MM)

global nn; global wm; global hm; global h; global ne; global G; global w_speed;
global Le1; global Le2; global Le3; global Le4; global Le5; global Le6;

GC = sparse(zeros(nn,nn));

%% RIGHT
CeR =  [0    0    0    0;
    0    0.6667    0.3333    0;
    0    0.3333    0.6667    0;
    0    0    0    0];

for e = [wm/h:wm/h:ne]
    el1 = ismember(e,Le1);
    if el1 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(1)/w_speed(1)) * (h/2) * CeR(i,j)+ GC(mi,mj);
            end
        end
    end
    el2 = ismember(e,Le2);
    if el2 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(2)/w_speed(2)) * (h/2) * CeR(i,j)+ GC(mi,mj);
            end
        end
    end
    el3 = ismember(e,Le3);
    if el3 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(3)/w_speed(3)) * (h/2) * CeR(i,j)+ GC(mi,mj);
            end
        end
    end
    el4 = ismember(e,Le4);
    if el4 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(4)/w_speed(4)) * (h/2) * CeR(i,j)+ GC(mi,mj);
            end
        end
    end
    el5 = ismember(e,Le5);
    if el5 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(5)/w_speed(5)) * (h/2) * CeR(i,j)+ GC(mi,mj);
            end
        end
    end
    el6 = ismember(e,Le6);
    if el6 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(6)/w_speed(6)) * (h/2) * CeR(i,j)+ GC(mi,mj);
            end
        end
    end
end



%% LEFT
CeL =  [0.6667    0    0    0.3333;
    0    0    0    0;
    0    0    0    0;
    0.3333    0    0    0.6667];


for e = [1:wm/h:ne]
    el1 = ismember(e,Le1);
    if el1 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(1)/w_speed(1)) * (h/2) * CeL(i,j)+ GC(mi,mj);
            end
        end
    end
    el2 = ismember(e,Le2);
    if el2 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(2)/w_speed(2)) * (h/2) * CeL(i,j)+ GC(mi,mj);
            end
        end
    end
    el3 = ismember(e,Le3);
    if el3 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(3)/w_speed(3)) * (h/2) * CeL(i,j)+ GC(mi,mj);
            end
        end
    end
    el4 = ismember(e,Le4);
    if el4 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(4)/w_speed(4)) * (h/2) * CeL(i,j)+ GC(mi,mj);
            end
        end
    end
    el5 = ismember(e,Le5);
    if el5 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(5)/w_speed(5)) * (h/2) * CeL(i,j)+ GC(mi,mj);
            end
        end
    end
    el6 = ismember(e,Le6);
    if el6 == 1
        for i = 1:4
            for j = 1:4
                mi = MM(e,i);
                mj = MM(e,j);
                GC(mi,mj) = (G(6)/w_speed(6)) * (h/2) * CeL(i,j)+ GC(mi,mj);
            end
        end
    end
end


%% BOTTOM
CeB =  [0    0    0    0;
    0    0    0    0;
    0    0    0.6667    0.3333;
    0    0    0.3333    0.6667];


for e = [1:wm/h]
    for i = 1:4
        for j = 1:4
            mi = MM(e,i);
            mj = MM(e,j);
            GC(mi,mj) = (G(1)/w_speed(1)) * (h/2) * CeB(i,j)+ GC(mi,mj);
        end
    end
end
