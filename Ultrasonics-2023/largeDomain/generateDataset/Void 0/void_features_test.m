function [] = void_features_test()
%% Void features:
global L; global H; global Lpml;
global N_void; global Eli_fea; global xRotated_total; global yRotated_total;

%%
N_void = 2;                                      % Creates number of voids of random number.
if N_void==0
    Eli_fea = [0,0,0,0,0];
else
    for p= 1: N_void
        loc= randi([1 2]);
        Y_loc= [0.15 0.35];
        
        a= 0.03+ (0.08-0.03).*rand(1,1);                                             % Major axis= m/2
        b= 0.01+ (0.015-0.01).*rand(1,1);                                            % Minor axis= n/2
        H_ellipse=  a+ ((L- a)-a).*rand(1,1) ;                                 % Random center x ellipse
        K_ellipse= Y_loc(loc); % 0.15;  %H/2 + ((H- a)-H/2).*rand(1,1) ;                                  % Random center y ellipse
        Angle= -5 + (5-(-5)).*rand(1,1);                                       % Random angle ellipse
        Eli_fea(p,:)=[a,b,H_ellipse,K_ellipse,Angle];       % Saving the features of ellipse
        
        
        q= linspace(0, 360,1000);
        transformMatrix = [cosd(Angle), sind(Angle); -sind(Angle), cosd(Angle)];
        xAligned =   a * sind(q);
        yAligned =   b * cosd(q);
        xyAligned = [xAligned; yAligned]';
        xyRotated = xyAligned * transformMatrix;
        xRotated = H_ellipse+ xyRotated(:, 1);
        yRotated = K_ellipse+ xyRotated(:, 2);
        xRotated_total(:,p)=xRotated;                      % Shape of each ellipse in x
        yRotated_total(:,p)=yRotated;                      % Shape of each ellipse in y
        
    end
end