clear all
close all
clc

%choose which algorithm to run
% alg = {BBOA, DCS, FS}

alg = 'BBOA';

%set nonmonotone memory size M (=1 for monotone)
M = 4;

%set the problem to be solved
pname     = 'davidon 2 (b)';
pdim      = 4; % problem has 4 variables
m         = 2; % problem has 2 general constraints
lbint     = zeros(pdim,1);
ubint     = 100.0*ones(pdim,1);
x_initial = 50.0*ones(pdim,1);
startp    = [25.0;5.0;-5.0;-1.0];
lb        = startp-10.0;
ub        = startp+10.0;
fhandle   = @(x)davidon2_b(x,lb,ub,lbint,ubint);
max_fun   = 5000;
outlev    = 0;

outline = ['Solving problem ' pname ' using algorithm '];
if M > 0
    outline = [outline 'NM-' alg ': '];
else
    outline = [outline ' M-' alg ': '];
end
fprintf('\n\n%35s',outline);

%call the optimizer
[x,f, stopfl, Dused] = box_DFL(alg, fhandle, m, M, x_initial, lbint, ubint, max_fun, outlev);
