function [Point_loading_x] = point_loading_2(Ricker_freq)

global Ntimestep; global dt;

%%
for i = 1:2*Ntimestep+1
    t = (i-1)*(dt/2);
    
    wc = 2*pi*Ricker_freq;
    tbar = 6*sqrt(6)/wc;
    
    if t <= tbar
        eta = wc*t-3*sqrt(6);
        y = -5000*((0.25*eta^2-0.5).*exp(-0.25*eta.^2)-13.*exp(-13.5))...
            /(0.5+13.*exp(-13.5));
    else
        y = 0;
    end
    Point_loading_x(i,1) = y;
end
end