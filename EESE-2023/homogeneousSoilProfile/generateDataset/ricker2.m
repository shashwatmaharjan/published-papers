function y = ricker2(fc,A,t,delay)

%% Ricker
t2 = t-delay;
wc = 2*pi*fc;
tbar = 6*sqrt(6)/wc;

if t2 <= tbar

    if t2 < 0
        y = 0;
    else
        eta = wc*t2-3*sqrt(6);
        y = -A*((0.25*eta^2-0.5).*exp(-0.25*eta.^2)-13.*exp(-13.5))...
            ./(0.5+13.*exp(-13.5));
    end

else
    y = 0;
end

end