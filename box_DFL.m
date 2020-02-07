%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% box_DFL
% Copyright (C) 2018 G.Liuzzi, S.Lucidi, F.Rinaldi:
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <https://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,f, stopfl, Dused] = box_DFL(alg, func_f, mm, M, x_initial_in, lb, ub, max_fun, outlev, varargin)
%
% Function box_DFL solves the integer problem:
%
%                   min      f(x)
%                   s.t. lb <= x <= ub
%                            g_i(x) <= 0,  i=1,...,mm
%                              x \in Z^n
%
% The user must provide func_f to evaluate the function f, lower and upper
% bounds (lb and ub vectors)
%
% Inputs:
%
%         alg      : which algorithm to use. One of BBOA, DCS, FS
%
%         func_f   : handle of the function to be minimized, i.e. f(x).
%
%         mm       : number of constraints. If problem is unconstrained,
%                    then mm = 0
%
%         M        : dimension of the memory for non-monotone search (must be >= 1) 
%
%         x_initial: the initial point to start the optimizer.
%
%         lb, ub   : lower and upper bounds.
%
%         max_fun  : maximum number of allowed function evaluations
%
%         outlev   : 0 for no output, >0 otherwise
%
% Outputs:
%
%          x       : best point
% 
%          f       : o.f. value related to x
%
%          stopfl  : 1 if code stopped for alpha< threshold
%                    2 if code stopped for nf >= max_fun
%                   99 if initial point does not satisfy the bound
%                         constraints
%          Dused   : set of directions used in the last iteration
%
% Functions called: either funct or functpen (provided by the optimizer),
%                   nm_discrete_linesearch (provided by the optimizer),
%                   nm_discrete_search (provided by the optimizer),
%                   generate_dirs (provided by optimizer).
%
% Written by G. Liuzzi, S. Lucidi, F. Rinaldi, 2016.
%

format compact;
format long;

x_initial = round(x_initial_in);

ubint = round(ub);
lbint = round(lb);

if (sum(abs(ub-ubint)) ~= 0) || (sum(abs(lb-lbint)) ~= 0)
    fprintf('\n\nERROR: upper and/or lower bound on some variable is NOT integer.');
    fprintf('\n       Please correct and resubmit.\n\n');
    x       = x_initial;
    f       = Inf;
    stopfl  = 99;
    Dused   = [];
    return;
end

m = max(mm,0);
if M < 1
    M = 1;
end
if ~strcmp(alg,'DCS') && ~strcmp(alg,'FS')
    alg = 'BBOA';
end

iter        = 0;      % iteration counter
alpha_start = 1;      % the starting stepsize 
nf          = 0;      % number of function evaluations performed
cache_hits  = 0;      % number of times point found in cache
n           = length(x_initial); % dimension of the problem
stop        = 0;      % stopping condition
W           = NaN*ones(1,M); % memory of function value for NM linesearch
xW          = NaN*ones(n,M);

nnf=0;
xf=Inf*ones(max_fun,n+m+1);

Phalton     = haltonset(n);
ihalton     = 7;
eta         = 1.5;
allones     = 0;

rho         = 1.0;
eps         = 0.1*ones(m,1);

chk_feas = (x_initial >= lb) & (x_initial <= ub);
if min(chk_feas) == 0
    fprintf('\n\nInitial point does not satisfy the bound constraints!\n\n');
    x       = x_initial;
    f       = Inf;
    stopfl  = 99;
    Dused   = [];
    return;
end

% D           denotes the set of search direction (one per each column)
% alpha_tilde is a row vector of stepsizes (one per each direction in D)


x  = x_initial;
if(m > 0)
    f  = functpen(x);
else
    f  = funct(x);
end
nf = nf+1;

if(m > 0)
    g = xf(1,n+2:n+m+1);
    rho = max(1.0,sum(max(g,0)));
    eps(max(g,0)<1.0) = 1.e-3;
    f  = functpen(x);
end

W(1) = f;
xW(:,1) = x;
bestf = f;
bestx = x;

D = eye(n,n);
successes = zeros(1,n); % components of this vector counts number of successes each direction has had
                        % where "success" means that alpha > 0 is returned
                        % by LS along the direction itself

if strcmp(alg,'BBOA')
    alpha_tilde   = round((ub+lb)/2.0)';
else
    alpha_tilde   = ones(1,n);
end
old_maxalpha  = Inf;
ndir          = size(D,2);

%  | 12345 | 12345 | +1234567890123 | +1234567890123 | +1234567890123 |
%  |  iter |   nf  |       f        |       f_ref    |    max_alpha   |
if outlev > 0
    if(m > 0)
        print_format = ['| %5d | %5d | %5d | %+13.8e | %+13.8e | %+13.8e | %+13.8e | %5d |   \n'];
        fprintf('\n');
        fprintf('|  iter |    nf | cache |        fpen     |        f_ref    |         viol    |    max_alpha    |  ndir |\n');
    else
        print_format = ['| %5d | %5d | %5d | %+13.8e | %+13.8e | %+13.8e | %5d |   \n'];
        fprintf('\n');
        fprintf('|  iter |    nf | cache |        f        |        f_ref    |    max_alpha    |  ndir |\n');
    end
else
    fprintf('   fun.evals =      ');
end

while ~stop 
    iter = iter+1;
    y = x; fy = f;
    cache_hits = 0;
    for idir = 1:ndir
        d = D(:,idir);
        if(iter == 1)
            f_ref = W(1);
        else
            f_ref = max(W);
        end
        
        if strcmp(alg,'BBOA') || strcmp(alg,'DCS')
            [alpha, x_trial, f_trial] = nm_discrete_linesearch(y,d,alpha_tilde(idir),lb,ub,f_ref);
        else
            [alpha, x_trial, f_trial] =     nm_discrete_search(y,d,alpha_tilde(idir),lb,ub,f_ref);
        end
        if alpha <= 0
            d = -d;
            if strcmp(alg,'BBOA') || strcmp(alg,'DCS')
                [alpha, x_trial, f_trial] = nm_discrete_linesearch(y,d,alpha_tilde(idir),lb,ub,f_ref);
            else
                [alpha, x_trial, f_trial] =     nm_discrete_search(y,d,alpha_tilde(idir),lb,ub,f_ref);
            end
            if alpha > 0
                successes(idir) = successes(idir)+1;
                if allones >= 1
                    allones = 0;
                end
                D(:,idir) = d;
                y  = x_trial;
                fy = f_trial;

                alpha_tilde(idir) = alpha;
                W = circshift(W,[0,1]);
                xW= circshift(xW,[0,1]);
                W(1) = fy;
                xW(:,1) = y;
                if(fy < bestf)
                    bestf = fy;
                    bestx = y;
                end
            else
                alpha_tilde(idir) = max(1,floor(alpha_tilde(idir)/2));
            end
        else
            successes(idir) = successes(idir)+1;
            if allones >= 1
                allones = 0;
            end
            
            y  = x_trial;
            fy = f_trial;
            alpha_tilde(idir) = alpha;
            W = circshift(W,[0,1]);
            xW= circshift(xW,[0,1]);
            W(1) = fy;
            xW(:,1) = y;
            if(fy < bestf)
                bestf = fy;
                bestx = y;
            end
        end

        if m > 0
            if (allones >= 1)
               break
            end
        else
            if (allones > 1)
               break
            end
        end
    end
    
    if(m > 0)
            sxf=size(xf,1);
            diff=(xf(:,1:n)-repmat(y',sxf,1)).^2;
            [mn,ind]=min(sum(diff,2));
            g = xf(ind(1),n+2:n+m+1);
    end
    
    if (norm(y-x) <= 1.e-14) && (max(alpha_tilde) == 1) && (old_maxalpha == 1)

        allones=allones+1;

        if strcmp(alg,'DCS') || strcmp(alg,'FS')
            if(bestf < fy)
                y  = bestx;
                fy = bestf;
            else
                stopfl = 1;
                stop   = 1;
                Dused  = D;
            end
            iexit = 1;
        else
            iexit = 0;
        end
        
        while iexit == 0
            % enrich set D
            [D, successes, alpha_tilde, iexit] = generate_dirs(n,D,successes,alpha_tilde,eta,0);
            if iexit == 0
                eta = eta + 0.5;
            end
            if eta >= 0.5*(norm(ub - lb)./2)
                %stop execution
                if(bestf < fy)
                    y  = bestx;
                    fy = bestf;
                else
                    stopfl = 1;
                    stop   = 1;
                    Dused  = D;
                end
                iexit = 1;
            end
        end
        
        ndir    = size(D,2);

        if(m > 0)
            %check on the penalty parameters eps
            ind_change = (max(g,0)>rho);
            nchg = size(ind_change,1);
            eps_changed = 0;
            for i = 1:nchg
                if(eps(ind_change(i)) > 1.e-10)
                    eps(ind_change(i)) = eps(ind_change(i))/2.0;
                    allones = 0;
                    eps_changed = 1;
                end
            end
            if(eps_changed == 1)
                sxf=size(xf,1);
                diff=(xf(:,1:n)-repmat(x',sxf,1)).^2;
                [mn,ind]=min(sum(diff,2));
                if (mn<=10^-16)
                     fval = xf(ind(1),n+1);
                     gval = xf(ind(1),n+2:end);
                     gval = gval';
                     m = size(gval,2);
                     f = fval + sum(max(gval,0)./eps);
                     cache_hits = cache_hits+1;
                else
                     if m > 0
                         f = functpen(x);
                     else
                         f = funct(x);
                     end
                     nf = nf + 1;
                end
                W           = NaN*ones(1,4);
                xW          = NaN*ones(n,4);
                W(1)        = f;
                xW(:,1)     = x;
                bestf       = f;
                bestx       = x;
            end
        end
    end

    if m > 0
        rho = max(1.e-8,rho*0.5);
    else
        rho = rho/2.0;
    end
    
    if (norm(y-x) <= 1.e-14) && (ndir >=5000)
        stopfl = 1;
        stop   = 1;
        Dused  = D;
        x      = bestx;
        f      = bestf;
    end
    
    x = y; f = fy;

    old_maxalpha = max(alpha_tilde);
    
    if outlev > 0
        if(m > 0)
            fprintf(print_format, iter, nf, cache_hits, f, f_ref, sum(max(g,0)), max(alpha_tilde),ndir);
        else
            fprintf(print_format, iter, nf, cache_hits, f, f_ref, max(alpha_tilde),ndir);
        end
    else
        fprintf('\b\b\b\b\b%5d',nf);
    end
    
    if nf >= max_fun
        stopfl = 2;
        stop   = 1;
        Dused  = D;
        x      = bestx;
        f      = bestf;
    end
    
end

if outlev == 0
    fprintf('\n');
end

function floc = funct(xint)
    floc = feval(func_f,xint);

    nnf = nnf+1;
    nx=size(xint,1);
    xf(nnf,1:nx)=xint';
    xf(nnf,nx+1)=floc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF CODE funct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function fpen = functpen(xint)
    [floc, gloc] = feval(func_f,xint);

    nnf = nnf+1;
    nx=size(xint,1);
    ng=size(gloc,1);
    fpen = floc + sum(max(gloc,0)./eps);
    xf(nnf,1:nx)=xint';
    xf(nnf,nx+1)=floc;
    xf(nnf,nx+2:nx+ng+1) = gloc';
end

function [alpha, x, f] = nm_discrete_linesearch(y,d,alpha_tilde,lb,ub, f_ref)
    %
    % Function nm_discrete_linesearch
    %
    % Purpose:
    %
    % This function performs a nonmonotone discrete linesearch
    % along a given direction d (d \in Z^n)
    %
    % Inputs:
    %
    % y            : starting point for the linesearch
    %
    % d            : search direction
    %
    % alpha_tilde  : starting stepsize
    %
    % lb, ub       : lower and upper bounds
    %
    % f_ref        : reference o.f. value 
    %
    % Output:
    %
    %
    % alpha        : 1) alpha > 0 if linesearch finds a point guaranteeing 
    %                simple decrease: f(y+alpha d)<f_ref;
    %                2) alpha = 0 failure
    %
    % x            : best point found in the linesearch 
    %
    % f            : o.f. value related to x 
    %

    % calculate dimension of the problem
    n = length(d);

    % initialize vector alpha_max
    alpha_max = Inf * ones(n,1);

    % caluclate max alpha
    indices = ( d > 0 );

    alpha_max(indices)=( ub(indices) - y(indices) )./ d(indices);


    indices = ( d < 0 );

    alpha_max(indices)=( lb(indices) - y(indices) )./ d(indices);

    %compute starting alpha
    alpha_bar  = floor( min(alpha_max) );
    alpha_init = min(alpha_tilde, alpha_bar);

    %Build first point for starting linesearch
    if (alpha_init > 0)
        y_trial = y + alpha_init * d;
        
        sxf=size(xf,1);
        diff=(xf(:,1:n)-repmat(y_trial',sxf,1)).^2;
        [mn,ind]=min(sum(diff,2));
        %diff
        %keyboard
        if (mn<=10^-16)
            fval = xf(ind(1),n+1);
            gval = xf(ind(1),n+2:end);
            gval = gval';
            if(m > 0)
                f_trial = fval + sum(max(gval,0)./eps);
            else
                f_trial = fval;
            end
            cache_hits = cache_hits+1;
        else
            if m > 0
                f_trial = functpen(y_trial);
            else
                f_trial = funct(y_trial);
            end
            nf = nf + 1;
        end
        
    else
        f_trial = Inf;
    end

    % cicle for updating alpha
    if (alpha_init > 0) && (f_trial<f_ref)

        % initialize alpha and best point
        alpha=alpha_init;
        x=y_trial;
        f=f_trial;

        %calculate trial point
        if(alpha < alpha_bar)
            y_trial = y + min(alpha_bar,2* alpha) * d;
            sxf=size(xf,1);
            diff=(xf(:,1:n)-repmat(y_trial',sxf,1)).^2;
            [mn,ind]=min(sum(diff,2));
            %diff
            %keyboard
            if (mn<=10^-16)
                fval = xf(ind(1),n+1);
                gval = xf(ind(1),n+2:end);
                gval = gval';
                if(m > 0)
                    f_trial = fval + sum(max(gval,0)./eps);
                else
                    f_trial = fval;
                end
                cache_hits = cache_hits+1;
            else
                if m > 0
                    f_trial = functpen(y_trial);
                else
                    f_trial = funct(y_trial);
                end
                nf = nf + 1;
            end

        else
            f_trial = Inf;
        end

        % expansion step (increase stepsize)
        while (alpha<alpha_bar) && (f_trial < f_ref)

            % alpha calulation and best point updating
            alpha=min(alpha_bar, 2*alpha);

            % best point updating
            x=y_trial;
            f=f_trial;


            %next point to be tested
            if(alpha < alpha_bar)
                y_trial = y + min(alpha_bar, 2* alpha) * d;
                sxf=size(xf,1);
                diff=(xf(:,1:n)-repmat(y_trial',sxf,1)).^2;
                [mn,ind]=min(sum(diff,2));
                %diff
                %keyboard
                if (mn<=10^-16)
                    fval = xf(ind(1),n+1);
                    gval = xf(ind(1),n+2:end);
                    gval = gval';
                    if(m > 0)
                        f_trial = fval + sum(max(gval,0)./eps);
                    else
                        f_trial = fval;
                    end
                    cache_hits = cache_hits+1;
                else
                    if m > 0
                        f_trial = functpen(y_trial);
                    else
                        f_trial = funct(y_trial);
                    end
                    nf = nf + 1;
                end
            else
                f_trial = Inf;
            end

        end
    else
        alpha = 0;
        x = y;
        f = Inf;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF CODE nm_discrete_linesearch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [alpha, x, f] = nm_discrete_search(y,d,alpha_tilde,lb,ub, f_ref)
    %
    % Function nm_discrete_search
    %
    % Purpose:
    %
    % This function performs a nonmonotone discrete linesearch
    % along a given direction d (d \in Z^n)
    %
    % Inputs:
    %
    % y            : starting point for the linesearch
    %
    % d            : search direction
    %
    % alpha_tilde  : starting stepsize
    %
    % lb, ub       : lower and upper bounds
    %
    % f_ref        : reference o.f. value 
    %
    % Output:
    %
    %
    % alpha        : 1) alpha > 0 if linesearch finds a point guaranteeing 
    %                simple decrease: f(y+alpha d)<f_ref;
    %                2) alpha = 0 failure
    %
    % x            : best point found in the linesearch 
    %
    % f            : o.f. value related to x 
    %

    % calculate dimension of the problem
    n = length(d);
    
    % initialize vector alpha_max
    alpha_max = Inf * ones(n,1);
    
    % caluclate max alpha
    indices = ( d > 0 );
    
    alpha_max(indices)=( ub(indices) - y(indices) )./ d(indices);
    
    
    indices = ( d < 0 );
    
    alpha_max(indices)=( lb(indices) - y(indices) )./ d(indices);
    
    %compute starting alpha
    alpha_bar  = floor( min(alpha_max) );
    alpha_init = min(alpha_tilde, alpha_bar);
    
    %Build first point for starting linesearch
    if (alpha_init > 0)
        y_trial = y + alpha_init * d;
        
        sxf=size(xf,1);
        diff=(xf(:,1:n)-repmat(y_trial',sxf,1)).^2;
        [mn,ind]=min(sum(diff,2));
        %diff
        %keyboard
        if (mn<=10^-16)
            fval = xf(ind(1),n+1);
            gval = xf(ind(1),n+2:end);
            gval = gval';
            if(m > 0)
                f_trial = fval + sum(max(gval,0)./eps);
            else
                f_trial = fval;
            end
            cache_hits = cache_hits+1;
        else
            if m > 0
                f_trial = functpen(y_trial);
            else
                f_trial = funct(y_trial);
            end
            nf = nf + 1;
        end
        
        if (f_trial<f_ref)
            x=y_trial;
            f=f_trial;
            alpha=alpha_init;
        else
            f_trial = Inf;
            alpha = 0;
            x = y;
            f = Inf;
        end
        
    else
        f_trial = Inf;
        alpha = 0;
        x = y;
        f = Inf;
    end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF CODE nm_discrete_search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function flag = prime_vector(d)
    n = length(d);
    flag = 0;
    if(n==1)
        flag = 1;
        return;
    end
    temp = gcd(abs(d(1)),abs(d(2)));
    if(n==2)
        flag = (temp == 1);
        return;
    end
    for i = 3:n
        temp = gcd(temp,abs(d(i)));
        if temp == 1
            flag = 1;
            return;
        end
    end
    if temp ~= 1
        flag = 0;
        return
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF CODE prime_vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function [Dout, succout, alpha, iexit] = generate_dirs(n,D,succ,alpha_tilde,eta,betaLS)
    %
    % Function generate_dirs
    %
    % Purpose:
    %
    % This function generate new integer directions which are added to set D 
    %
    % Inputs:
    %
    % n            : dimension of the problem
    %
    % D            : matrix of current directions (one per each column) 
    %
    % alpha_tilde  : array of stepsizes along direction already in D
    %
    % Output:
    %
    % Dout         : [new_direction D] 
    %
    % succout      : [0 succ] 
    %
    % alpha        : array of stepsizes along the directions in Dout
    %                alpha = [new_step_sizes alpha_tilde]
    %

    % d = rand(n,1);
    % d = d./norm(d);
    % 
    % Q = [null(d') d];
    % 
    % Dout = [Q D];
    % alpha = [ones(1,n) alpha_tilde];

    mD = size(D,2);

    for j = 1:1000
        %keyboard
        v = 2.0*Phalton(ihalton,:)' - ones(n,1);
        ihalton = ihalton+1;
        v = eta*v./norm(v);

        if (norm(v)<1.e-16)
           break
        end

        %d = abs(round(v)); good if H=norm(d)^2*eye(n,n) - 2*d*d' used
        d = round(v);



        %now check whether d is a prime vector
        if prime_vector(d)
            trovato = 0;
            %check whether d is already in D
            DIFF1=D-repmat(d,1,mD);
            DIFF2=D+repmat(d,1,mD);
            if( min ( sum(abs(DIFF1),1)) == 0 ) || ( min ( sum(abs(DIFF2),1)) == 0 ) 
                trovato = 1;    
            end

            if ~trovato
                H       = [d]; %norm(d)^2*eye(n,n) - 2*d*d';
                Dout    = [H D];
                succout = [0 succ];
                alpha   = [max(betaLS, max(alpha_tilde)) alpha_tilde]; %ones(1,n)]; 
                iexit   = 1;
                return
            end
        end
    end

    Dout    = D;
    succout = succ;
    alpha   = alpha_tilde;
    iexit   = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF CODE generate_dirs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF CODE box_DFL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end