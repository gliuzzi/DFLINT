function y = kowalik(xint,lb,ub,lbint,ubint)
% 
% Nonsmooth KOWALIK-OSBORNE function
% Matlab Code by G.Liuzzi (October 4, 2017).
% The number of variables n should be adjusted below.
% The default value of n = 2.
% 

x = lb + ((ub - lb)./(ubint - lbint)).*xint;

z = [0.1957 4.0;
     0.1947 2.0;
     0.1735 1.0;
     0.1600 0.5;
     0.0844 0.25;
     0.0627 0.1670;
     0.0456 0.1250;
     0.0342 0.1;
     0.0323 0.0833;
     0.0235 0.0714;
     0.0246 0.0625];
u = z(:,2);

f = (x(1)*(u.^2 + x(2)*u))./(u.^2 + x(3)*u + x(4)) - z(:,1);

y = max(abs(f));

end
