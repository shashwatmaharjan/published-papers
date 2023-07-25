% load('Signal_data_RK.mat')
load('Signal_ani_void_1_1.mat')
% load('Crack_data_void_1.mat')

for i= 1: length(um_history)
     for sensors= 1:41

         A= um_history{i,1}(sensors, :);
         check_nan = isnan(A);
         for j= 1: length(check_nan)
            if check_nan(j) == 1
             store (i,1) = 1;
            elseif check_nan(j) == 0
                store(i,1) = 0;
            end
         end

     end
     
end

count= 1;
for e= 1: length(um_history)
    if store(e,1)== 1
        delete(count,1)= e;
        count= count + 1;
    end
end

for i= length(delete):-1:1
    um_history(delete(i,1),:)= [];
%     Crack_data_yes_no_1(delete(i,1),:)= [];
end

% save('Crack_data_void_1.mat','Crack_data_yes_no_1')
% save('Signal_data_ricker.mat','um_history')