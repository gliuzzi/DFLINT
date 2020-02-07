clear all
close all
clc

%choose which algorithm to run
% alg = {BBOA, DCS, FS}

alg = 'BBOA';

%set nonmonotone memory size M (=1 for monotone)
M = 4;

%set the problem to be solved
pname     = 'kowalik-osborne';
pdim      = 4; % problem has 4 variables
m         = 0; % problem has 0 general constraints
lbint     = zeros(pdim,1);
ubint     = 100.0*ones(pdim,1);
x_initial = (ubint+lbint)/2;
startp    = [0.25; 0.39; 0.415; 0.39];
lb        = startp-10.0;
ub        = startp+10.0;
fhandle   = @(x)kowalik(x,lb,ub,lbint,ubint);
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
