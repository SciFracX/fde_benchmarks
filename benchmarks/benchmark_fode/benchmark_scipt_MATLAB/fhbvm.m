function [t,y,stats,err] = fhbvm( fun, y0, T, M, txton )
%
%   [t,y[,stats[,err]]] = fhbvm( fun, y0, T[, M] )
%
%   Solves the Caputo fractional IVP:
%
%    y^(alpha) = f(t,y),   0<=t<=T,   
%
%    y^(i)(0)  = y_i \in R^m, i = 0,...,[alpha]-1,
%
%   with [alpha] = ceil(alpha).  N.B.: alpha>0 and not an integer.
%
%    Input:
%    ======
%    - y0    : matrix [alpha] x m with the initial conditions. 
%              In more detail,
%
%                     y0(i,:) = y_{i-1}^T,  i = 1,...,[alpha].
%
%    - T     : final integration time.
%
%    - fun   : a string containing a function identifier, 
%              or a function_handle, such that:
%
%              alpha   = fun()
%
%              f(t,y)  = fun( t, y )   (N.B.: MUST WORK  IN  VECTOR MODE,
%                                             i.e., if [ell,m] = size(y),
%                                             and  length(t) == ell, then
%                                             fun(t,y)   must  return  an  
%                                             ell x m  matrix)
%
%              f'(t,y) = fun( t, y, 1 ) (i.e., computes the Jacobian of f).
%
%
%    - M     : an optional parameter such that h = T/M is the
%              stepsize of a uniform mesh, if possible, or the
%              largest stepsize (approximately) of a graded mesh.
%              If M is not given in input, a default value
%              between 2 and 10 is guessed. In general, M 
%              must be >1 and should be as small as possible.
%
%    - txton : an optional flag to have text output and warnings
%              (this affects the execution time).
%
%    Output:
%    =======
%    - (t,y) : computed solution;
%
%    - stats : vector containing the following times:
%              - pre-processing time, 
%              - problem solving time, 
%              - pre-processing time for error estimation, 
%              - problem solving time for error estimation;
%  
%    - err   : (optional) if required in output, the global error is 
%              estimated on a doubled mesh (this is quite costly).
%
%                                                          Rel. 2025-05-04.
%
%    References:
%    ===========
%    [1] L.Brugnano, K.Burrage, P.Burrage, F.Iavernaro. 
%                 A spectrally accurate step-by-step method for the 
%                 numerical solution of fractional differential equations. 
%                 J. Sci. Comput. 99 (2024) 48.
%                 https://doi.org/10.1007/s10915-024-02517-1
%    [2] L.Brugnano, G.Gurioli, F.Iavernaro. 
%                 Numerical solution of FDE-IVPs by using Fractional HBVMs: 
%                 the fhbvm code. Numerical Algorithms 99 (2025) 463-489.
%                 https://doi.org/10.1007/s11075-024-01884-y
%    [3] L.Brugnano, G.Gurioli, F.Iavernaro. 
%                 Solving FDE-IVPs by using Fractional HBVMs: some 
%                 experiments with the fhbvm code. 
%                 J. Comput. Methods Sci. Eng. 25(1) (2025) 1030-1038.
%                 https://doi.org/10.1177/14727978251321328 
%                 (preprint on arXiv.2407.11460 [math.NA])
%    [4] Website: https://people.dimai.unifi.it/brugnano/fhbvm/
%
if nargout==4
   Nmax = 2.5e3;  %% maximum number of mesh points with error estimation %%
else
   Nmax = 5e3;    %%%% maximum number of mesh points %%%%
end   
if nargin<4 
   M  = min( ceil(T/0.5), 10 ); 
   M  = max( M, 3 );   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2025-04-10
elseif M<=1 || M~=fix(M)
   error('fhbvm: an integer M > 1 is needed')
end
if T<=0
   error('fhbvm: T>0 is needed')
end
if isempty(fun) || ( ~ischar(fun) && ~isa(fun,'function_handle') ) 
   error('fhbvm: fun is wrong') 
end    
alfa      = feval(fun);
if alfa<=0 || alfa==round(alfa) 
   error('fhbvm: a non integer and positive alpha is required')
elseif alfa>size(y0,1)
   error('fhbvm: alpha cannot exceed the number of initial conditions')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   initialize execution time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stats     = [0 0 0 0];
err       = [];
data      = fhbvmdata();
data.alfa = alfa;
data.fun  = fun;
data.T    = T;
data.y0   = y0;
data.n    = M;
if nargin==5 
   data.txt = 1;
else
   data.txt = 0; warning off
end   
data.k    = 22;         k    = data.k;
data.s    = 20;         s    = data.s; 
[x,beta,PO,Is,~,Xs,abd] = jacobidat( k, s, alfa );
x         = x(:);
x1        = [x; 1];     % the abscissa 1 is also needed to advance the step
gam       = findgam( Xs, alfa ); 
CC        = gam*inv(Xs);
fattb     = min( norm(PO,1)*norm(Is,1), norm(PO,inf)*norm(Is,inf) );
data.x    = x;
data.x1   = x1;
data.beta = beta;
data.PO   = PO;
data.Is   = Is;
data.abd  = abd;
data.fatt = fattb;
data.gam  = gam;
data.CC   = CC;
data      = findh0( data );
if data.h0>T, error('fhbvm: h0 cannot be larger than T'), end
h0        = data.h0;
n         = data.n;   
r         = data.r;  
Ps        = data.Ps;     
glc       = data.glc;    
glb       = data.glb;  
if n>Nmax 
   warning(['number of mesh-points too high; '...
            'the problem could be not solved!']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% For a graded mesh,
%
%  Jj(:,nu*k+i) = Jj(0:s-1,(r^nu-1)/(r-1)+c_ir^nu), nu = 1,...,n-1,
%
% or, for a uniform mesh,
%
%  Jj(:,nu*k+i) = Jj(0:s-1,nu+c_i), nu = 1,...,n-1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Jj     = zeros(s,(k+1)*(n-1));
ind    = 0;
for nu = 1:n-1
    for i = 1:k+1
        ind = ind+1;
        if r>1 % graded mesh
           M = (r^nu-1)/(r-1)+x1(i)*r^nu;
        else   % uniform mesh
           M = nu+x1(i);
        end   
        Jj(:,ind) = Jjalfa1( abd, x, beta, alfa, M, Ps, glc, glb );
    end    
end
data.Jj  = Jj;
tcor     = toc;
stats(1) = tcor;  
[t,y,~,~,~,data] = fhbvmdata( data );         % data.NaN is set
stats(2) = toc-tcor; 
tcor     = toc;  
if data.NaN       %%%%%%%%%%%%%% an error message is issued %%%% 2025-04-08
   fprintf(' NaN occurred at t = %g \n', t(end) ); 
elseif nargout==4 %%%%%%%%%%%%%% error estimation
   r2       = sqrt(r);
   n2       = 2*n;
   if r2==1
      h02   = h0/2;
   else    
      h02   = h0*(r2-1)/(r-1);   
   end    
   data.r   = r2; 
   data.n   = n2;
   data.h0  = h02;
   data.err = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   r = r2;  n = n2;   % because of the doubled mesh %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Jj     = zeros(s,(k+1)*(n-1));
   ind    = 0;
   for nu = 1:n-1
       for i = 1:k+1
           ind = ind+1;
           if r>1 % graded mesh
              M = (r^nu-1)/(r-1)+x1(i)*r^nu;
           else   % uniform mesh
              M = nu+x1(i);
           end   
           Jj(:,ind) = Jjalfa1( abd, x, beta, alfa, M, Ps, glc, glb );
       end    
   end
   data.Jj  = Jj;
   stats(3) = toc-tcor; 
   tcor     = toc;  
   [t2,y2,~,~,~,data] = fhbvmdata( data );       %%%%%%%%%%%%%%% 2025-05-04
   if data.NaN                                    
      disp('NaN in error estimation'), err = []; % NaN: no error estimation
   else                                          % OK: error estimation       
      err      = zeros(size(y2));
      dt       = norm(t-t2(1:2:end),inf);
      if dt>10*n2*eps*(1+T)
         warning(['fhbvmdata: error estimation maybe inaccurate,' ...
                  ' dt = %g \n'], dt )
      end    
      err(1:2:end,:)   = abs(y-y2(1:2:end,:));
      err(2:2:end-1,:) = ( err(1:2:end-2,:)+err(3:2:end,:) )/2;      % mean 
      t        = t2;
      y        = y2;                % the doubled-mesh solution is returned
   end   
   stats(4) = toc-tcor;  
end    
return

function [t,y,etim,it,blend,data] = fhbvmdata( data )  
%
%  [t,y,etim,it,blend,data] = fhbvmdata( data );
%
%     or
%
%  data = fhbvmdata;   (initialize a void structure for later use).
%
%  Solves the system of fractional ODEs, for alfa\in(0,1), 
%
%         D^alfa y(t) = f(t,y), t \in [0,T], 
%            y^(i)(0) = y_0^i \in R^m, i = 0,...,[alfa]-1,         
%
%  with D^alfa the Caputo fractional derivative, and [alfa] the ceil of 
%  alfa.  The problem is solved by using the FHBVM(22,20) method,  with 
%  blended iteration when necessary.
%  If required, the global error for is estimated on a doubled mesh.
%  (the default is not, since it is costly).
%
%
%  Input:
%---------------------------------------
%
%   data - structure containing all data needed by the method (see below, 
%          except data.NaN):
%
%            data.txt  - flag for display (default==1);
%            data.err  - flag for error estimation (default==0);
%            data.fun  - function identifier, string or function_handle, 
%                        mandatory; fun - function definig the problem:
%
%                        alfa    = fun 
%
%                        f(t,y)  = fun(t,y)   (has to work in vector mode);
%
%                        f'(t,y) = fun(t,y,1)  (vector mode not needed);
%
%            data.k    - parameter k, set at k = 22;
%            data.s    - parameter s, set at s = 20 of FHBVM(k,s);
%            data.alfa - order alfa of the derivative. If not set, it is
%                        evaluated by fun (0<alfa).
%                     >> The output value cannot be modified for subsequent 
%                        calls;
%            data.y0   - initial condition. If not set, it is evaluated
%                        by fun;
%            data.T    - final time T. If not set, it is evaluated by fun.
%                     >> The output value cannot be modified for subsequent 
%                        calls;
%            data.n    - total number of mesh-point, besides t_0=0; 
%            data.h0   - initial stepsize h0; 
%            data.r    - parameter for the graded mesh:
%                        t_0 = 0, 
%                        t_i = t_(i-1)+h_(i-1) == t_(i-1)+r^(i-1)*h0, 
%                                                          i=1,...,n;
%            data.x    - Jacobi abscissae;
%            data.x1   - x1 := [x; 1];
%            data.beta - Jacobi weights;
%            data.PO   - matrix P_s^T Omega of FHBVM; 
%            data.Is   - matrix Is^alfa of FHBVM;
%            data.abd  - coefficients for the evaluation of Jacobi
%                        polynomials;
%            data.fatt - |PO|*|Is|;
%            data.gam  - parameter for the blended iteration;
%            data.CC   - matrix CC := gam*( PO*Is )^(-1);
%            data.Jj   - auxiliary integrals for the memory terms. 
%                        If not given, they are evaluated;
%            data.Ps   - Jacobi polynomials evaluated at the GL abscissae;
%            data.glc  - GL abscissae of enough high order;
%            data.glb  - GL weights; 
%            data.NaN  - flag set to .TRUE. if a NaN has occurred
%                        (only needed in output): it is the only field set
%                        by fhbvmdata;
%
%--------------------------------------------------------------------------
%
%  Output:
%---------------------
%    (t,y) - solution;
%
%     etim - vector of length 2 with the elapsed times:
%            etim(1) = time for the memory terms of the local problems;
%            etim(2) = time for solving the local problems;
%
%       it - vector of length 2 with the required iterations:
%            it(1) = total number of blended iterations for solving the
%                    local problems;
%            it(2) = total number of fixed-point iterations for solving the
%                    local problems;
%
%    blend - number of timesteps where the blended iteration has been used
%            (the remaining ones have been solved by using the fixed-point
%             iteration);
%
%     data - structure containing all auxiliary fractional integrals and
%            matrices. Only data.NaN is set on output;
%
%     fin  - norms of the fundamental matrix at the mesh points (optional).
%
%--------------------------------------------------------------------------
%                                                          Rel. 2025-04-07.
%--------------------------------------------------------------------------
if nargin==0
   data.txt  = []; 
   data.fun  = [];
   data.k    = 22;   % fixed at this moment
   data.s    = 20;   % fixed at this moment
   data.alfa = [];
   data.y0   = [];
   data.T    = [];
   data.n    = [];
   data.h0   = [];
   data.r    = [];
   data.x    = [];
   data.x1   = [];
   data.beta = [];
   data.PO   = [];
   data.Is   = [];
   data.abd  = [];
   data.Ps   = [];
   data.glc  = [];
   data.glb  = [];
   data.fatt = [];
   data.gam  = []; 
   data.CC   = []; 
   data.Jj   = []; 
   data.NaN  = [];
   t         = data;
   return
end   
%--------------------------------------------------------------------------
etim     = [0 0]; 
txt      = data.txt;       % flag for selected output display, if set to 1.  
fun      = data.fun;
k        = data.k;
s        = data.s;
alfa     = data.alfa;
T        = data.T;
n        = data.n;         % N.B.: n has to be >1, in case of a graded mesh
y0       = data.y0;
h0       = data.h0;
r        = data.r;
x        = data.x;
x1       = data.x1;
PO       = data.PO;
Is       = data.Is;
fattb    = data.fatt;
gam      = data.gam;
CC       = data.CC;
Jj       = data.Jj; 
m        = size(y0,2);
h        = h0*cumprod( [1; r*ones(n-1,1)] );     %  stepsizes
t        = [0; cumsum(h)];
if abs(T-t(end))>max(100,2*n)*(1+T)*eps
   disp('----------------------------------------------')
   disp('fhbvmdata: if r==1, T==n*h0;                  ') 
   disp('           if r>1,  T==h0*(r^n-1)/(r-1).      ')
   disp('Currently, n, r, h0, and T are not compatible.')
   disp('----------------------------------------------')
   error('fhbvmdata: exit')
else                                             % correct the mesh
   cT = T/t(end);
   h  = h*cT;
   t  = t*cT; 
end    
y        = zeros(n+1,m);
y(1,:)   = y0(1,:);
fatt     = 1/gamma(alfa);
fatt1    = fatt/alfa;
inds     = 0;
indk     = 0;
s1       = 1:s;
halfa    = h.^alfa;
it       = zeros(1,2);
%
%  gamJ(:,(nu-1)*s:nu*s) = Jacobi-Fourier coefficients at step nu
%                          for the given problem
gamJ     = zeros(m,n*s); 
Im       = eye(m);
%-------------------------------------------------------------------------%
if txt==1
   linea = '-------------------------------------------------------------';
   disp(linea)
   func = fun;
   if ~ischar(func), func = func2str(func); end
   disp( ['Problem ' func ' solved with:'] );
   fprintf( 'FBHVM(%2i,%2i), alfa = %g \n', k, s, alfa );
   if r==1
      disp( 'on a uniform mesh' ); 
   else
      disp( 'on a graded mesh' ); 
   end   
end
%-------------------------------------------------------------------------%
flag_NaN = 0;
small    = 0.6;  % bound for using the blended iteration instead of fix-pt.
blend    = 0;    % number of blended steps
for i = 1:n       
    etcor           = toc;  
    J0              = feval( fun, t(i), y(i,:), 1 ); 
    L               = min( norm(J0,1), norm(J0,inf) );
    blendflag       = L*fattb*halfa(i) > small;
    if blendflag
       teta         = inv( Im -halfa(i)*gam*J0.' );    % transpose Jacobian
       blend        = blend + 1;
    else
       teta         = [];
    end   
    tcor            = t(i) + h(i)*x;
    tcor1           = [tcor; t(i+1)];
    fi0             = evalfi( y0, tcor1, x1, Jj(:,1:indk), gamJ(:,1:inds));
    etim(1)         = etim(1) + toc - etcor;     % memory term problem time
    etcor           = toc;
    [yi,gami,iti,flag_NaN] = step( fun, tcor, fi0, halfa(i)*Is, PO, ...
                                             fatt1*halfa(i), CC, teta, L );
    etim(2)         = etim(2) + toc - etcor;     % solution problem time      
    if blendflag 
       it(1)        = it(1) + iti;
    else   
       it(2)        = it(2) + iti;
    end   
    if flag_NaN
       warning('fhbvmdata: NaN at step %i, exit \n', i)
       t = t(1:i);    y = y(1:i,:);
       break
    end    
    y(i+1,:)        = yi;
    gamJ(:,inds+s1) = ( fatt*halfa(i) )*gami.';
    inds            = inds + s;
    indk            = indk + k+1;
end    
data.NaN = flag_NaN;
%-------------------------------------------------------------------------%
if txt==1
   disp(linea)
   fprintf( ' %i blended steps and %i fixed-point steps\n', ...
                                              blend, i-blend ); 
   fprintf( '       total iterations solver: \n' );
   fprintf( '                      blended = %i \n', it(1) );
   fprintf( '                  fixed-point = %i \n', it(2) );
   disp(linea)
   fprintf( 'r = %g, n = %i, interval = [0,%g], h0 = %g \n', ...
                                             r, i, t(end), h0 );
   if flag_NaN
      fprintf(' NaN occurred at step %i \n', i );
      if errflag~=0
         disp(' no error estimation done');
      end    
   end   
   disp(linea)
   fprintf( '                  memory time = %g \n', etim(1) );
   fprintf( '                  solver time = %g \n', etim(2) );
      disp( '                  ------------------------    ' );
   fprintf( '                   TOTAL TIME = %g \n', sum(etim) );
   disp(linea)
end
%-------------------------------------------------------------------------%
return

function [data1,tim2] = findh0( data )
%
%   Computes the parameters h0, n, and r of a graded mesh, besides other 
%   informations re-used by the main program.
%
%                                                          Rel. 2025-04-10.
m         = length(data.y0);
T         = data.T;
N         = data.n;
abd       = data.abd;
x         = data.x;
x1        = data.x1;
beta      = data.beta;
alfa      = data.alfa; 
k         = data.k;
s         = data.s;
h         = T/N;
tol       = m*300*eps;
Nmax      = 5;     % maximum value of N, to allow a uniform mesh for ell=2.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data1     = data;
data1.Jj  = [];       % this has to be void
data1.n   = 1;        % only 1 step is done
data1.T   = h;
data1.h0  = h;
data1.r   = 1;
data1.txt = 0;                               
imax      = 10;
fatN      = 1;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2025-04-09 
for i = 1:imax
    [~,y,~,~,~,data1] = fhbvmdata( data1 );   
    if data1.NaN==0
       break           % successfull step
    else
       fatN     = fatN*2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 2025-04-09
       h        = h/2; % the stepsize is halved
       data1.T  = h;
       data1.h0 = h;   % disp([i imax])
       if i==imax, warning('findh0: initial step may be wrong'), end
    end
end    
yend      = y(end,:);                        
ellmax    = 33;    % corresponds to decrease h by a factor < 1e-19                       
data1.r   = 3;  % value of r for obtaining a reduction factor 4 per iterate
data1.n   = 2;  % number of sub-intervals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tcor      = toc;
Ps        = [];  % Jacobi polynomials aveluated at the GL abscissae
c         = [];  %     Gauss-Legendre abscissae, and
b         = [];  %     Gauss-Legendre weights of suitable order
Jj        = zeros(s,(k+1));
r         = data1.r;
for i = 1:k+1
    M                = 1 + x1(i)*r;
    [Jj(:,i),Ps,c,b] = Jjalfa1( abd, x, beta, alfa, M, Ps, c, b );
end    
tim2      = toc-tcor;  %  pre-processing time
data1.Jj  = Jj;
data1.Ps  = Ps;   %
data1.glc = c;    %  for later use in main
data1.glb = b;    % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
for ell   = 1:ellmax
    data1.T           = h;
    data1.h0          = h/4;   % current estimate of the initial stepsize
    [~,y,~,~,~,data1] = fhbvmdata( data1 );   
    err               = norm( ( y(end,:)-yend )./( 1 + abs(yend) ), inf );  
    if err<=tol && data1.NaN==0, break, end  %%%%%%%%%%%%%%%%%%%%% 25-04-10
    yend              = y(2,:);
    h                 = h/4;   
end
if ell >= ellmax 
   if data1.NaN
      error('findh0: the initial timestep cannot be coputed')  % 2025-04-10
   else   
      warning('findh0: h0 possibly not appropriate'); h = h*4; % 2024-07-11
   end   
end
if     ell<=1  % a uniform mesh with N points is used
   r      = 1; 
   n      = N;
elseif ell<=2 && N<=Nmax % a uniform mesh with 4*N points is used
   r      = 1; 
   n      = 4*N;
elseif N==1
   error('findh0: N cannot be equal to 1 for a graded mesh')
else       % a graded mesh is used
   N      = N*fatN; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 25-04-09
   r      = ( N - 4^(1-ell) )/( N - 1 );      % initial guess for r      
   n      = ceil( 1 + log( 4^(ell-1) )/log(r) );     % number of timesteps  
  [r,~]   = findr( 1, 4^(1-ell)/N, n, r );    % final value of r
end
data1.h0  = h;
data1.Jj  = [];       % need to be recomputed for the final mesh
data1.n   = n;        % final number of timesteps
data1.T   = data.T;   % final time restored
data1.r   = r;        % updated value of r
data1.txt = data.txt; % restored original value
return 

function [r,flagr] = findr( T, h0, n, r0 )
%
% Compute the proper value of r for the graded mesh.
%
% Input:
%   T  = width of the integration interval,
%   h0 = initial stepsize,
%   n  = number of timesteps,
%   r0 = initial guess for r.
%
% Output:
%   r     = value of the parameter r,
%   flagr = true, if the final parameter is wrong.
%                                                          Rel. 2024-02-23.
if h0<=0, error('findr: h0 is wrong'), end
if T<=h0, error('findr: T has to be larger than h0'), end
if r0<=1, error('findr: initial guess for r <= 1'), end
Th0   = T/h0;
r     = r0;      
r0    = 1;
ii    = 0;
imax  = 1000;
while r0~=r && ii<imax
    ii = ii+1; r0 = r; r = ( 1 + (r-1)*Th0 )^(1/n);     
    if ii >= imax, disp([ii r r-r0]), break, end
end
flagr = r<1 || r0<1 || abs(r-r0)>eps*(r+r0);
if flagr
   warning('findr: parameter r for the graded mesh may be wrong')
end   
return

function gam = findgam( Xs, alfa )
%
%  Computes the optimal value for the blended iteration
%
%                                                          Rel. 2024-03-14.
lam = sort(eig(Xs));
gam = abs(lam);
s   = length(lam);
ro  = zeros(s,1);
for i = 1:s
    ro(i) = max( abs( lam-gam(i) ).^2./( 2*gam(i)*gam ) );
end
[rom,ind] = min(ro);
gam       = gam(ind);
if rom>1 && alfa<=1
   warning('findgam: blended iteration not A-convergent')
end    
return

function fi0 = evalfi( y0, t1, x1, Jj, gamJ )
%
%  Evaluates the memory term for the local FDE.
%
%                                                          Rel. 2025-03-07.
[s,nk1] = size(Jj);
[nu,m]  = size(y0);
k1      = length(x1);
fi0     = 0;  
fatj    = 1;
t1j     = ones(k1,1);  % t^j/j! , j=0,1,...
for j = 1:nu
    fi0  = fi0 + ( (t1j/fatj) )*y0(j,:);
    if j<nu, t1j = t1j.*t1; fatj = fatj*j; end
end    
n       = nk1/k1;
if n==0, return, end  %%%%%%%%%%%%%%%%%%%%%%%%%%>>>>>>  initial step, quit.

fi1  = zeros(k1,m);
indg = (n-1)*s;
indJ = 0;
s1   = 1:s;
kk1  = 1:k1;
for nu = 1:n
    fi1  = fi1 + ( gamJ(:,indg+s1)*Jj(:,indJ+kk1) )';
    indg = indg - s;
    indJ = indJ + k1;
end
fi0 = fi0 + fi1;
return

function [y,gam,it,fNaN] = step( fun, t, fi0, Is, PO, ghalfa, CC, teta, L )
%
%   Integration step for the local problem.
%                                                          Rel. 2024-02-08.
%
[k,m]  = size(fi0);
k      = k-1;
Y0     = fi0(1:k,:);
tol0   = 2*eps*(1+L);
tol    = 100*tol0;
itmax  = 1000*(m+1)*(k+1);
y      = Y0;
scal   = 1./( 1 + abs(Y0) );
err    = Inf;
gam    = 0;
noblen = isempty(teta); 
for it = 1:itmax
    yold = y;
    f    = feval(fun,t,y);
    if noblen
       gam  = PO*f;
    else
       eta  = gam - PO*f;
       eta1 = CC*eta;
       gam  = gam -( eta1 +(eta-eta1)*teta )*teta;
    end   
    y    = Y0 + Is*gam;
    err0 = err;
    err  = norm( (y-yold).*scal, inf );    
    if isnan(err) 
        warning('step: NaN'), break
    elseif err<=tol0 
        break 
    elseif err>err0 && err0<=tol
        break
    end  
end
fNaN = any(any(isnan(y)));  % a scalar flag is needed to check NaN
if ~fNaN 
   if err>tol*2 
      warning('step: step non converging: error = %g \n',err) 
   end
   y = fi0(k+1,:) + ghalfa*gam(1,:);
end   
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         AUXILIARY FUNCTIONS                             %   
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,beta,PO,Is,Ps,Xs,abd] = jacobidat( k, s, alfa )
%
%   [x,beta,PO,Is,Ps,Xs,abd] = jacobidat(k,s,alfa)
%
%   Gauss-Jacobi weights, abscissae, and matrices for FHBVM(k,s)
%
%  INPUT:
%       - k,s  : parameters for HBVMF(k,s), 1<=s<=k;
%       - alfa : index of the fractional equation ( alfa > 0 );
%       - prec : .true. for vpa, otherwise double.
%  OUTPUT:
%       - x    : k Gauss-Jacobi abscissae in [0,1], relative to the 
%                weight function
%                               w(x) = alfa*(1-x)^alfa;
%       - beta : corresponding quadrature weights; % sum(beta)=1 %
%       - PO   : PO = P_s^T * Omega, with P_s(i,j) = P_{j-1}(x_i)
%                                    and Omega = diag(beta);
%       - Is   : matrix I_s(i,j) = I^alfa P_{j-1}(x_i);
%       - Ps   : matrix P_s(i,j) = P_{j-1}(x_i);
%       - Xs   : = PO*Is;
%       - abd  : coefficients [a b d] for the recurrence:
%               P_{j+1}(c)=(a(j)*c-b(j))*P_j(c)-d(j)*P_{j-1}(c), j=0,..s-2,
%               P_0(c)==1, P_{-1}(c)==0.
%
%                                                          Rel. 2024-07-17.
%
if (s>k) || (s<1) || (alfa<=0) || (alfa==round(alfa)) 
    error('jacobidat: wrong parameters') 
end
ab   = r_jacobi01(k,alfa-1,0);
xw   = gauss(k,ab);
x    = xw(:,1);
beta = xw(:,2);
if nargout<=2, return, end
Ps   = zeros(k,k+1);
Ps(:,1) = 1;
Ps(:,2) = (x-ab(1,1)).*Ps(:,1);
for j = 2:k
    Ps(:,j+1) = (x-ab(j,1)).*Ps(:,j) - ab(j,2)*Ps(:,j-1);
end
err = norm(Ps(:,k+1),inf);  
if err>k*eps   % check that the Jacobi abscissae are accurate %
    warning('jacobidat: are abscissae accurate?'),  keyboard 
end
Ps = Ps(:,1:s);
% beta is normalized such that sum(beta)=1
beta = alfa*beta; % beta = beta/sum(beta);  %%%%% validated 2023.08.21
% polynomials are scaled in order to be orthonormal
scal = 1./sqrt(diag(Ps.'*diag(beta)*Ps));
Ps   = Ps*diag(scal);
PO   = Ps.'*diag(beta);

Is   = zeros(k,s);
for i = 1:k    %, disp([i k])
    Is(i,:) = makeint( ab(:,1), ab(:,2), alfa, x(i), s, x, beta );  
end    
Is   = Is*diag(scal);
Xs   = PO*Is;
%%%%%%%%%%%%%%%%%%%%%%%%%%%  validated 2023.08.21
abd  = zeros(s-1,3);
dx   = diag(x);
for j = 1:s-1
    aj = 1/( PO(j+1,:)*diag(x)*Ps(:,j) );
    bj = aj*PO(j,:)*dx*Ps(:,j);
    if j==1
        dj = 0;
    else
        dj = aj/abd(j-1,1);
    end    
    abd(j,:) = [aj bj dj];
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%
return

function Ialfa = makeint( b, c, alfa, x, s, ci, wi )
%
%    Ialfa = makeint( b, c, alfa, x, s, ci, wi )   
%
%    Compute I^alfa P_j(x), j = 0,...,s-1, 
%
%    I^alfa being the Caputo integral, 0<alfa,
%
%    where:
%           P_0(x) = 1
%           P_1(x) = (x-b_1)*P_0(x)
%           P_j(x) = (x-b_j)*P_{j-1}(x) -c_j*P_{j-2}(x),  j=2,...,s-1.
%
%    Input:
%        (b,c) - coefficients for building the monic polynomials (c(1)=0);
%      alfa, x - order and ending point of the integral.
%
%                                                          Rel. 2024-07-17.
%
k  = length(ci);
Ps = zeros(k,s);
xx = x*ci;
Ps(:,1) = 1;
if s>1, Ps(:,2) = xx - b(1); end   
for j = 2:s-1
    Ps(:,j+1) = ( xx - b(j) ).*Ps(:,j) - c(j)*Ps(:,j-1);
end    
Ialfa = ( x^alfa/gamma(alfa+1) )*( wi.'*Ps )';
return

function [Ialfa,Ps,x,beta] = Jjalfa1( abd, cJ, wJ, alfa, M, Ps, x, beta )
%
%   [Ialfa,Ps,x,beta] = Jjalfa1( abd, cJ, wJ, alfa, M, Ps, x, beta )  
%
% Computes the integrals:
% 
%       J_j^alfa(M) = int_0^1(M-c)^(alfa-1)P_j(c)dc, j = 0...s-1,
%
% with  {Pj(c)} satisfying the recurrence:
%
%       P_0(c) = 1
%       P_1(c) = (a(1)*c-b(1))
%       P_j(c) = (a(j)*c-b(j))*P_{j-1}(c)-d(j)*P_{j-2}(c), j = 2...s-1,
%
% and abd = [a b d].  
%
% If x>=Mmax, a Gauss-Legendre formula of order 2*k is used, with nodes and 
% weights (x,beta) (k=30).    Conversely, a recursion formula is used.
%
% Ps contains the Jacobi polynomials (with parameters (1-alfa,0)) evaluated
% at the Gauss-Legendre points x.
% 
% Ps, x, beta are computed, if void in input.
%
%                                                          Rel. 2024-01-31.
%
Mmax = 1.1;  
k    = 30;      % order of the Gauss-Legendre quadrature = 2*k
if isempty(Ps)          
   a        = abd(:,1);
   b        = abd(:,2);
   d        = abd(:,3);
   s        = size(abd,1)+1;
   [x,beta] = GL(k);  % Gauss-Legendre abscissae and weights
   x        = x(:);
   beta     = beta(:);
   Ps       = zeros(k,s);
   Ps(:,1)  = 1;
   if s>1 
      Ps(:,2)  = a(1)*x-b(1); 
      for j = 2:s-1
          Ps(:,j+1) = (a(j)*x-b(j)).*Ps(:,j) -d(j)*Ps(:,j-1);
      end
   end   
end
if M>=Mmax
   Ialfa = Ps.'*( beta.*( (M-x).^(alfa-1) ) );  % Gauss-Legendre quadrature
else
   Ialfa = Jjalfa( abd, cJ, wJ, alfa, M );     
end                                                  
return

function Ialfa = Jjalfa( abd, ci, wi, alfa, x )
%
% Ialfa = Jjalfa( abd, ci, wi, alfa, x )
%
% Computes the integrals:
%
%       J_j^alfa(x) = int_0^1(x-c)^(alfa-1)P_j(c)dc, j = 0...s-1,
%
% by using the Gauss-Jacobi quadrature (ci,wi), with the polynomials
% {Pj(c)} satisfying the recurrence:
%
%       P_0(c) = 1
%       P_1(c) = (a(1)*c-b(1))
%       P_j(c) = (a(j)*c-b(j))*P_{j-1}(c)-d(j)*P_{j-2}(c), j = 2...s-1.
%
% N.B.: abd  = [a b d], and x \in [1,1.1), for better performance.
%
%
%                                                          Rel. 2024-01-20.
%
s       = size(abd,1)+1; 
k       = length(ci);
a       = abd(:,1);
b       = abd(:,2);
d       = abd(:,3);
Ps      = zeros(k,s);
Ps(:,1) = 1;
if s>1
   xx      = x*ci;
   Ps(:,2) = a(1)*xx - b(1);
end   
for j = 2:s-1
    Ps(:,j+1) = ( a(j)*xx - b(j) ).*Ps(:,j) -d(j)*Ps(:,j-1);
end
Ialfa   = x^alfa*( wi.'*Ps );
Ps(:,1) = 1;
if s>1
   xx      = 1 + (x-1)*ci;
   Ps(:,2) = a(1)*xx - b(1);
end   
for j = 2:s-1
    Ps(:,j+1) = ( a(j)*xx - b(j) ).*Ps(:,j) -d(j)*Ps(:,j-1);
end
Ialfa = Ialfa - (x-1)^alfa*( wi.'*Ps );
Ialfa = Ialfa(:)/alfa;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   Functions from: 
%     OPQ: A MATLAB SUITE OF PROGRAMS FOR GENERATING ORTHOGONAL POLYNOMIALS
%     AND RELATED QUADRATURE RULES.
%   Companion Matlab codes of:
%     W.Gautschi, Orthogonal Polynomials Computation and Approximation,
%                 Oxford University Press, 2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ab = r_jacobi01(N,a,b)
%R_JACOBI01 Recurrence coefficients for monic shifted Jacobi polynomials.
%   AB=R_JACOBI01(N,A,B) generates the Nx2 array AB of the first
%   N recurrence coefficients for monic shifted Jacobi polynomials 
%   orthogonal on [0,1] relative to the weight function 
%   w(x)=(1-x)^A*x^B. The call AB=R_JACOBI01(N,A) is the same as
%   as ab=R_JACOBI01(N,A,A) and AB=R_JACOBI01(N) the same as
%   AB=R_JACOBI01(N,0,0).
%
%   Based on: 
%     OPQ: A MATLAB SUITE OF PROGRAMS FOR GENERATING ORTHOGONAL POLYNOMIALS
%     AND RELATED QUADRATURE RULES.
%   Companion Matlab codes of:
%     W.Gautschi, Orthogonal Polynomials Computation and Approximation,
%                 Oxford University Press, 2004.
if nargin<2, a = 0; end;  if nargin<3, b = a; end
if((N<=0)||(a<=-1)||(b<=-1))
    error('r_jacobi01: parameter(s) out of range') 
end
cd = r_jacobi(N,a,b);
n  = 1:N;  ab(n,1) = (1+cd(n,1))./2;  
ab(1,2) = cd(1,2)/2^(a+b+1);
n = 2:N;  ab(n,2) = cd(n,2)./4;
return

function xw = gauss(N,ab)
%GAUSS Gauss quadrature rule.
%   GAUSS(N,AB) generates the Nx2 array XW of Gauss quadrature
%   nodes and weights for a given weight function W. The nodes,
%   in increasing order, are placed into the first column of XW,
%   and the corresponding weights into the second column. The
%   weight function W is specified by the Nx2 input array AB
%   of recurrence coefficients for the polynomials orthogonal
%   with respect to the weight function W.
%
%   Based on: 
%     OPQ: A MATLAB SUITE OF PROGRAMS FOR GENERATING ORTHOGONAL POLYNOMIALS
%     AND RELATED QUADRATURE RULES.
%   Companion Matlab codes of:
%     W.Gautschi, Orthogonal Polynomials Computation and Approximation,
%                 Oxford University Press, 2004.
N0 = size(ab,1); 
if N0<N, error('gauss: input array ab too short'), end
J  = zeros(N);
for n = 1:N, J(n,n) = ab(n,1); end
for n = 2:N
  J(n,n-1) = sqrt(ab(n,2));
  J(n-1,n) = J(n,n-1);
end
[V,D] = eig(J);
[D,I] = sort(diag(D));
V     = V(:,I);
xw    = [D ab(1,2)*V(1,:)'.^2];
return

function ab = r_jacobi(N,a,b)
%R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
%   AB=R_JACOBI(N,A,B) generates the Nx2 array AB of the first
%   N recurrence coefficients for the monic Jacobi polynomials 
%   orthogonal on [-1,1] relative to the weight function
%   w(x)=(1-x)^A*(1+x)^B. The call AB=R_JACOBI(N,A) is the same
%   as AB=R_JACOBI(N,A,A) and AB=R_JACOBI(N) the same as
%   AB=R_JACOBI(N,0,0).
%
%   Supplied by Dirk Laurie, 6-22-1998; edited by Walter
%   Gautschi, 4-4-2002.
%   Edited by Walter Leopardi 10-22-2006.
%
%   Based on: 
%     OPQ: A MATLAB SUITE OF PROGRAMS FOR GENERATING ORTHOGONAL POLYNOMIALS
%     AND RELATED QUADRATURE RULES.
%   Companion Matlab codes of:
%     W.Gautschi, Orthogonal Polynomials Computation and Approximation,
%                 Oxford University Press, 2004.
if nargin<2, a = 0; end  
if nargin<3, b = a; end
if((N<=0)||(a<=-1)||(b<=-1)) 
    error('r_jacobi: parameter(s) out of range') 
end
nu = (b-a)/(a+b+2);
if a+b+2 > 128
  mu = exp((a+b+1)*log(2)+((gammaln(a+1)+gammaln(b+1))-gammaln(a+b+2)));
else
  mu = 2^(a+b+1)*((gamma(a+1)*gamma(b+1))/gamma(a+b+2));
end
if N==1 
   ab = [nu mu]; 
else
   N  = N-1; n = 1:N; nab = 2*n+a+b;
   A  = [nu (b^2-a^2)*ones(1,N)./(nab.*(nab+2))];
   n  = 2:N; nab = nab(n);
   B1 = 4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
   B  = 4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
   ab = [A.' [mu; B1; B.']];
end   
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Other functions for computing the Gauss-Legendre weights and abscissae  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,beta] = GL(k)
%
%  [x,beta] = GL(k)
%
%  Gauss-Legendre abscissae and weights of order 2k.
%
%                                                          Rel. 2023-08-21.
%
[x,beta,~] = lgwt01(k); 
[x,ind]    = sort(x); 
beta       = beta(ind); 
return

function [x,w,y] = lgwt01(N)
%  
% This script computes the Legendre-Gauss nodes (x) and weights 
% (w) on the interval [0,1] with truncation order N.
% Also the nodes (y) on [-1,1] are computed.
% Adapted from the file written by Greg von Winckel-25/02/2004.
% 
%                                                          Rel. 2024-07-17.
%
a  = 0; b = 1;  % normalization on the interval [0,1] %
N  = N-1;
N1 = N+1; N2 = N+2;
xu = linspace(-1,1,N1)';
% Initial guess %
y  = cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
L  = zeros(N1,N2);   % Gauss-Legendre Vandermonde Matrix %
%%%Lp = zeros(N1,N2);   % Derivative of the GL-VM %    non serve 2024-07-17
% Compute the zeros of the N+1 Legendre Polynomial using %
% the recursion relation and the Newton-Raphson method   %
y0 = 2;                    % Iterate until new points are uniformly %
while max(abs(y-y0))>eps   %  within epsilon of old points          %    
    L(:,1) = 1;  L(:,2) = y; 
    %%%Lp(:,1) = 0; Lp(:,2) = 1;  non servono 2024-07-17
    for k = 2:N1
        L(:,k+1) = ( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    Lp = (N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    y0 = y;  y = y0-L(:,N2)./Lp;
end
x = (a*(1-y)+b*(1+y))/2;  % tranform from [-1,1] to [a,b]==[0,1] %    
w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2; % Compute the weights %
w = ( w + w(end:-1:1) )/2;              % and symmetrize them %
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF