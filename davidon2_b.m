function [y, g] = davidon2_b(xint,lb,ub,lbint,ubint)
%----------------------------
% Function Davidon 2
%----------------------------
    x = lb + ((ub - lb)./(ubint - lbint)).*xint;

    for i = 1:21
        t(i) = 0.25 + 0.75*(i-1)/20.0;
        f(i) = x(4) - (x(1)*t(i)^2 + x(2)*t(i) + x(3))^2 - sqrt(t(i));
    end
    y = max(abs(f));
    J = 1:(length(x)-2);
    g = (3-2*x(J+1)).*x(J+1) - x(J) - 2*x(J+2) + 2.5;
end
