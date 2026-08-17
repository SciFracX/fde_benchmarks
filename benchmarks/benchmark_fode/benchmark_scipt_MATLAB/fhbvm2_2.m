function [t,y,etim] = fhbvm2_2( fun, y0, T, N, n, nu, txton )
%--------------------------------------------------------------------------
%
%   [t,y,etim] = fhbvm2_2( fun, y0, T, N, n, nu )
%
%   Solves the initial value problem for fractional differential equations 
%   of Caputo type:
%
%    y1^(alpha1) = f1(t,y),      0<=t<=T,   
%    y2^(alpha2) = f2(t,y),
%
%    y^(i)(0)    = y_i \in R^m,  i = 0,...,ell-1,   (ell>=1),
%
%    with
%
%         y = [y1 y2] \in IR^m,  y1 \in IR^m1,  y2 \in IR^(m-m1), 
%
%   and   ell-1 < alpha1, alpha2 < ell  (possibly equal).
%
%   REMARK: if alpha2 is void or <=0, alpha1==alpha2 is assumed;
%           in such a case, m1 is not needed (it is set to m).
%
%   It uses a FHBVM(k,s) method using:
%
%     - the Jacobi abscissae, if alfa1==alfa2, with k = s = 22;
%     - the Jacobi-Pineiro abscissae, otherwise, with k = 30 and s = 22.
%
%    Input:
%    ======
%    - y0    : matrix  ell x m  with the initial conditions. 
%              In more detail,
%
%                     y0(i+1,:) = (y_i).',  i = 0,...,ell-1. 
%
%    - T     : final integration time.
%
%    - fun   : a string containing a function identifier, 
%              or a function_handle, such that:
%
%              [alphas,m1] = fun()
%
%              with alphas = alpha1
%                or alphas = [alpha1 alpha2]
%                 
%                   REMARK:   if alpha2 is void or <=0, then alpha1==alpha2
%                   -------   is assumed and, if so, m1 is not used.
%
%              f(t,y)  = fun( t, y ) (REMARK: MUST WORK  IN  VECTOR MODE,
%                                             i.e., if [k,m] = size(y),
%                                             and  length(t) == k, then
%                                             fun(t,y)   must return  a  
%                                             k x m  matrix)
%
%              f'(t,y) = fun( t, y, 1 ) (i.e., computes the Jacobian of f
%                                              at (t,y)).
%
%
%    - N     : h = T/N is the stepsize of the uniform mesh.
%
%    - n     : graded mesh in [0,n*h].
%
%    - nu    : number of graded mesh points in the interval  [0,n*h]  with 
%                        r = min( 2, n/(n-1) ). 
%              Consequently,  h1 = n*T*(r-1) / (N*(r^nu-1)) is the initial
%              step. This allows for coping for possible non smooth vector
%              fields, at t=0. 
%      REM: if  n==nu,  a uniform mesh is used; if  n==N, a graded mesh is
%           used; if  1<n<N  and  nu>n,  a mixed mesh is used.
%
%    Output:
%    =======
%    - (t,y) : computed solution;
%
%    - etim  : elapsed time.
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
%    [4] L.Brugnano, G.Gurioli, F.Iavernaro, M.Vikerpuur. 
%                 Analysis  and  implementation of  collocation methods for
%                 fractional differential equations.
%                 J. Sci. Comput. 104 (2025) 92.
%                 https://doi.org/10.1007/s10915-025-03006-9
%    [5] L.Brugnano, G.Gurioli, F.Iavernaro, M.Vikerpuur.
%                 A multi-order extension of Fractional HBVMs (FHBVMs). 
%                 J. Sci. Comput. 107 (2026) 94.
%                 https://doi.org/10.1007/s10915-026-03315-7
%    [6] Website: https://people.dimai.unifi.it/brugnano/fhbvm/
%
%
%                                                          Rel. 2026-05-25.
%
%--------------------------------------------------------------------------
warning on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% flag for the final output
txt       = nargin==7; 
data.txt  = txt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1 = n;    %%%%%  REMARK: hereafter  n1<-n  and  n<-N. N and n were  %%%%%%
n  = N;    %%%%%     used to be consistent with the notation in [4]. %%%%%%             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Checking the input   
if nargin<6, error('insufficient number of parameters'), end              %
if T<=0, error('wrong final time T'), end                                 %
if n~=fix(n) || n1~=fix(n1) || nu~=fix(nu), error('N,n,nu integer'), end  %
if n<n1, error('N>=n required'), end                                      % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% important check 
if nu<=n1, nu = 1; n1 = 1; if txt, disp('uniform mesh set'), end, end     %
uno = nu==n1;  % flag: 1, for a pure uniform mesh, 0 otherwise %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialize structure with code data  
data.fun   = fun;
s          = 22;  %%%%%%%%%%%%%%%%%%%%%%%%%% order of the method in [s,2*s] 
data.s     = s;
alfa       = feval(fun);  %%%%%%%%%%%% alpha(s)  
if isempty(alfa)
    error('void alpha(s) returned from the input function')
elseif length(alfa)>1
    alfa1  = alfa(1);
    alfa2  = alfa(2); 
    [~,m1] = feval(fun);  %%%%%%%%%%%% m1
else 
    alfa1  = alfa(1);
    alfa2  = []; 
end
if isempty(alfa2) || alfa2<=0 
    alfa2  = alfa1; 
end  
data.alfa1 = alfa1;
data.alfa2 = alfa2;
data.y0    = y0;
[ell,m]    = size(y0);
data.m     = m;
data.ell   = ell;
if min(alfa1,alfa2)<=ell-1 || max(alfa1,alfa2)>=ell  
    error( ' %i < alpha1, alpha2 < %i required', ell-1, ell )
elseif alfa1==alfa2
    flagalfa = 0; k = s; m1 = m;   % collocation at Jacobi abscissae
else
    flagalfa = 1; k = 30;          % Jacobi-Pineiro abscissae and weights
    if m1<=0 || m1>=m, error( 'm1 is wrong' ), end
end    
data.k     = k;
data.m1    = m1;     % length(m1)
m2         = m-m1;   
data.m2    = m2;     % length(y2)
%%%%%%%%%%%%%%%%%%%%%% 
data.T     = T;       %%% width of the integration interval
data.n     = n;
data.n1    = n1;
if n1==1, r = 2; else, r = n1/(n1-1); end
data.r     = r;       % parameter for the graded mesh
h          = T/n;     % timestep in the uniform mesh
data.h     = h;   
numin      = ceil(log(11)/log(r));  %%%%%%% guarantees hnu <= 1.1*h %%%%%%%
nup        = numin>nu && nu>n1;
if nup
    nu     = numin; 
end    
data.nu    = nu;                    %%% number of points in the graded mesh
data.h0    = n1*h*(r-1)/(r^nu-1);   %%% initial timestep in the graded mesh
h0         = data.h0;
hnu        = r^(nu-1)*h0;           %%% last timestep in the graded mesh
data.x     = [];
data.x1    = [];
%%%%%%%%%%%%%%%%%%%%%%%
data.beta1 = [];
data.PO1   = [];
data.Is1   = [];
data.abd1  = [];
data.Ps1   = [];
data.beta2 = [];
data.PO2   = [];
data.Is2   = [];
data.abd2  = [];
data.Ps2   = [];
%%%%%%%%%%%%%%%%%%%%%%%
data.glc   = [];
data.glb   = [];
data.fatt  = [];
data.Jj1   = []; 
data.Jjc1  = [];
data.Jj2   = []; 
data.Jjc2  = [];
data.NaN   = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialize execution time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the abscissae, weights, etc.
if flagalfa  
  %  
  % alfa1~=alfa2  (k=30, s=22):
  % ---------------------------
  %               Jacobi-Pineiro abscissae and corresponding weights 
  %               for quadratures of accuracy [1.5*k]=45 >= 44=2*s
  %
  [x,beta1,beta2,ier]      = myJacobiPineiro( k, alfa1, alfa2 );              
  if ier~=0 
     error( 'Jacobi-Pineiro abscissae and weights not properly computed' ) 
  end
  [~,~,PO1,Is1,~,Xs1,abd1] = jacobidat1( k, s, alfa1, x, beta1 );             
  [~,~,PO2,Is2,~,Xs2,abd2] = jacobidat1( k, s, alfa2, x, beta2 );
  data.gam = [];  % parameters for the blended 
  data.CC  = [];  % iteration not needed
else
  %   
  % alfa1==alfa2  (k=s=22):
  % -----------------------
  %               collocation at the Gauss-Jacobi abscissae with
  %               quadrature of accuracy 2*k=44 == 2*s
  %
  [x,beta1,PO1,Is1,~,Xs1,abd1] = jacobidat1( k, s, alfa1 );
  gam      = findgam( Xs1, alfa1 ); % parameters
  CC       = gam*inv(Xs1);          % for the
  data.gam = gam;                   % blended
  data.CC  = CC;                    % iteration
  beta2    = []; %
  PO2      = []; % 
  Is2      = []; % not needed when alfa1==alfa2
  Xs2      = []; %
  abd2     = []; %
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x          = x(:);
x1         = [x; 1];         % the abscissa 1 is needed to advance the step
if flagalfa
   X12     = PO1*Is2; % needed for the simplified
   X21     = PO2*Is1; % Newton iteration
   fattb   = min( norm(X12,1)+norm(X21,1)+norm(Xs1,1)+norm(Xs2,1), ...
               norm(X12,inf)+norm(X21,inf)+norm(Xs1,inf)+norm(Xs2,inf) );
else
   X12     = []; % not needed for the
   X21     = []; % blended  iteration
   fattb   = min( norm(PO1,1)*norm(Is1,1), norm(PO1,inf)*norm(Is1,inf) );
end   
data.fatt  = fattb;   % parameter to decide using the fixed-point iteration
data.x     = x;
data.x1    = x1;
%%%%%%%%%%%%%%%%%%%%
data.X11   = Xs1;
data.X12   = X12;
data.X22   = Xs2;
data.X21   = X21;
data.beta1 = beta1;
data.PO1   = PO1;
data.Is1   = Is1;
data.abd1  = abd1;
data.beta2 = beta2;
data.PO2   = PO2;
data.Is2   = Is2;
data.abd2  = abd2;
%%%%%%%%%%%%%%%%%%%%
kGauss     = 30;  % order of the Gauss-Legendre quadrature used == 2*kGauss
[Ps1,Ps2,glc,glb] = makePs1( abd1, abd2, kGauss ); % Gauss-Legendre 
data.kGauss  = kGauss;                             % abscissae and weights,
data.Ps1     = Ps1;                                % and corresponding 
data.Ps2     = Ps2;                                % Jacobi matrices
data.glc     = glc;   %  Gauss-Legendre abscissae and weights 
data.glb     = glb;   %  of order 2*kGauss = 60
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the auxiliary integrals
%
% For the graded mesh,
%
%  Jj(:,j*k+i) = Jj(0:s-1,(r^j-1)/(r-1)+c_ir^j), j = 0,...,nu-2,
%
% and, for the uniform mesh,
%
%  Jjc(:,j*k+i) = Jjc(0:s-1,j+c_i), j = 0,...,n-n1-1,
%                                                              i = 1,...,k.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nu>1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% graded mesh
   Jj1 = zeros(s,(k+1)*(nu-1));
   if flagalfa
      Jj2 = zeros(s,(k+1)*(nu-1));
   else
      Jj2 = []; 
   end   
else
   Jj1 = [];
   Jj2 = [];
end   
ind   = 0;
for j = 1:nu-1  %% graded mesh
    for i = 1:k+1
        ind       = ind+1;
        M         = (r^j-1)/(r-1) + x1(i)*r^j;
        Jj1(:,ind) = Jjalfa1( M, abd1, x, beta1, alfa1, Ps1, glc, glb ); 
        if flagalfa
           Jj2(:,ind) = Jjalfa1( M, abd2, x, beta2, alfa2, Ps2, glc, glb );
        end   
    end    
end
data.Jj1  = Jj1;
data.Jj2  = Jj2;
if n-n1>1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uniform mesh mesh
   Jjc1 = zeros(s,(k+1)*(n-n1));
   if flagalfa
      Jjc2 = zeros(s,(k+1)*(n-n1));
   else
      Jjc2 = []; 
   end   
else
   Jjc1 = [];
   Jjc2 = [];
end   
ind   = 0;
for j = 1:n-n1    %% h1*(r^nu-1)/(r-1) == n1*h
    for i = 1:k+1
        ind         = ind+1;
        M           = j + x1(i);
        Jjc1(:,ind) = Jjalfa1( M, abd1, x, beta1, alfa1, Ps1, glc, glb );
        if flagalfa
          Jjc2(:,ind) = Jjalfa1( M, abd2, x, beta2, alfa2, Ps2, glc, glb );   
        end  
    end    
end
data.Jjc1 = Jjc1;
data.Jjc2 = Jjc2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the numerical solution
[t,y,it,blend,flag_NaN] = fhbvmdata( data );   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
etim = toc;               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% wall clock time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
if txt==1
   clc 
   linea = '-----------------------------------------------------------'; 
   disp(linea)
   fprintf( 'Solver: FHBVM(%i,%i) \n', k, s )
   if flagalfa
      fprintf( 'alpha1 = %g  \n', alfa1 )
      fprintf( 'alpha2 = %g  \n', alfa2 )
      fprintf( 'm1     = %i  \n', m1 )
      fprintf( 'm2     = %i  \n', m2 )
   else   
      fprintf( 'alpha  = %g \n', alfa1 )
      fprintf( 'm      = %i \n', m )      
   end
   disp(linea)
   if nup
      fprintf( 'nu updated to %i \n', nu )
   end    
   if nu>n1
      fprintf( 'initial graded mesh with %i points \n', nu )
      fprintf( '     initial timestep    = %10.2e \n', h0 )
      fprintf( '     final   timestep    = %10.2e \n', hnu )
      fprintf( '     grading parameter r = %g \n', r )
      fprintf( '\n' )
   end
   if n>n1
      fprintf( '%i uniform points with \n', n-n1+uno )
      fprintf( '             timestep  h = %10.2e \n', h )
   end   
   disp(linea)
   if flagalfa
      fprintf( ' %i Newton steps and %i fixed-point steps\n', ...
                                              blend, length(t)-blend-1 ); 
   else   
      fprintf( ' %i blended steps and %i fixed-point steps\n', ...
                                              blend, length(t)-blend-1 ); 
   end
   fprintf( '    total iterations solver: \n' );
   if flagalfa
      fprintf( '                      Newton = %i \n', it(1) );
   else
      fprintf( '                     blended = %i \n', it(1) );
   end       
   fprintf( '                 fixed-point = %i \n', it(2) );
   disp( ' ' )
   if flag_NaN
      fprintf( ' NaN occurred at t = %g \n', t(end) );
   else
      fprintf( ' execution regularly ended at t = %g \n', t(end) ) 
   end   
   disp(linea)
   fprintf( ' execution time = %g sec \n', etim );
   disp(linea)
elseif flag_NaN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% an error message is issued
      fprintf( ' NaN occurred at t = %g \n', t(end) );  
end
%-------------------------------------------------------------------------%
return

function [t,y,it,blend,flag_NaN] = fhbvmdata( data )  
%
%  [t,y,it,blend,flag_NaN] = fhbvmdata( data );
%
%  Solves the system of fractional ODEs, for alfa>0, 
%
%         D^alfa_i y_i(t) = f_i(t,y), t \in [0,T], 
%            y_i(0) = y_0^i \in R^m_i, i = 1,2,
%            y = [y1 y2],  f = [f1 f2].
%
%  with D^alfa the  Caputo  fractional derivative, and [alfa] the ceil 
%  of alfa. The problem is solved by using the FHBVM(k,s) method, with 
%  blended iteration when necessary.
%
%  Input:
%---------------------------------------
%
%   data - structure containing all data needed by the method:
%
%            data.fun  - function identifier, string or function_handle, 
%                        mandatory; fun - function definig the problem:
%
%                        f(t,y)  = fun(t,y)   (has to work in vector mode);
%
%                        f'(t,y) = fun(t,y,1)  (vector mode not needed);
%
%            data.k    - parameter k of FHBVM(k,s);
%            data.s    - parameter s of FHBVM(k,s);
%            data.y0   - initial condition(s);
%            data.m    - dimension of the problem;
%            data.ell  - number of initial conditions;
%            data.T    - final time T; 
%            data.n    - h = T/n gives the uniform mesh;
%            data.h    - h as above;
%            data.nu   - nu mesh points in the graded mesh on [0,n1*h];
%            data.n1   - h1 as above;
%            data.h0   - initial stepsize h0; 
%            data.r    - parameter for the graded mesh:
%                        t_0 = 0, 
%                        t_i = t_(i-1)+h_(i-1) == t_(i-1)+r^(i-1)*h0, 
%                                                         i=1,...,nu;
%            data.x    - Jacobi-Pineiro abscissae;
%            data.x1   - x1 := [x; 1];
%            data.glc  - Gauss-Legendre abscissae of enough high order;
%            data.glb  - Gauss-Legendre weights;
%            data.gam  - parameter of the blended iteration;
%            data.CC   - for the blended iteration;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% i = 1/2
%            data.alfa1/2 - order alfa of the derivatives;
%            data.m1/2    - dimension of the corresponding subsystems;
%            data.beta1/2 - Jacobi-Pineiro weights;
%            data.PO1/2   - matrix P_s^T Omega of FHBVM; 
%            data.Is1/2   - matrix Is^alfa of FHBVM;
%            data.abd1/2  - coefficients for the evaluation of Jacobi
%                           polynomials;
%            data.fatt1/2 - |PO1/2|*|Is1/2|;
%            data.Jj1/2   - auxiliary integrals for the memory terms on the 
%                           graded mesh; 
%            data.Jjc1/2  - auxiliary integrals for the memory terms on the
%                           uniform mesh;
%            data.Ps1/2 - Jacobi polynomials evaluated at the GL abscissae.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%            it(1) = total number of Newton iterations for solving the
%                    local problems;
%            it(2) = total number of fixed-point iterations for solving the
%                    local problems;
%
%    blend - number of timesteps where the Newton iteration has been used
%            (the remaining ones have been solved by using the fixed-point
%             iteration);
%
%   flag_NaN - flag set to 1 in case NaN occur.
%
%--------------------------------------------------------------------------
%                                                          Rel. 2025-10-03.
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fun      = data.fun;                 %
y0       = data.y0;                  %
k        = data.k;                   %
s        = data.s;                   %
alfa1    = data.alfa1;               %
alfa2    = data.alfa2;               %
T        = data.T;                   %
n        = data.n;                   %
n1       = data.n1;                  %
nu       = data.nu;                  %
h0       = data.h0;                  %
r        = data.r;                   %
h        = data.h;                   %
x        = data.x;                   %
x1       = data.x1;                  %
PO1      = data.PO1;                 %
PO2      = data.PO2;                 %
Is1      = data.Is1;                 %
Is2      = data.Is2;                 %
X11      = data.X11;                 %
X12      = data.X12;                 %
X21      = data.X21;                 %
X22      = data.X22;                 %
fattb    = data.fatt;                %
Jj1      = data.Jj1;                 %
Jj2      = data.Jj2;                 %
Jjc1     = data.Jjc1;                %
Jjc2     = data.Jjc2;                %
m1       = data.m1;                  %
m2       = data.m2;                  %
m        = data.m;                   %
ell      = data.ell;                 %
gam      = data.gam;                 %
CC       = data.CC;                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hh       = [h0*cumprod( [1; r*ones(nu-1,1)] ); h*ones(n-n1,1)]; % timesteps
t        = [0; cumsum(hh)];                                     % and  mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abs(T-t(end))>max(100,2*n)*(1+T)*eps      %
   disp( '-------------------------------' ) %
   disp( 'fhbvmdata: mesh not compatible.' ) %  
   disp( '-------------------------------' ) %  
   error( 'fhbvmdata: exit' )                %    mesh consistency check
else                                         %    and correction  
   cT = T/t(end);                            %
   h  = h*cT;                                %
   t  = t*cT;                                %
end                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y        = zeros(nu+n-n1+1,m);   % y(1,:)        initial condition 
y(1,:)   = y0(1,:);              % y(2:nu+1,:)   graded mesh 
                                 % y(nu+2:end,:) uniform mesh
%%%% fatt  --> fatt_1,  fatt_2                                    
fatt_1   = 1/gamma(alfa1);
fatt_2   = 1/gamma(alfa2);  
%%%% fatt1 --> fatt1_1, fatt1_2
fatt1_1  = fatt_1/alfa1;                                                   
fatt1_2  = fatt_2/alfa2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flagalfa = alfa1~=alfa2; % different alfas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% halfa --> halfa1, halfa2
halfa1   = hh.^alfa1;   % h_i^alfa1 on the graded mesh  
%%%% half  --> half1,  half2
half1    = h^alfa1;     % h_i^alfa1 on the uniform mesh                    
if flagalfa
  halfa2 = hh.^alfa2;   % h_i^alfa2 on the graded mesh  
  half2  = h^alfa2;     % h_i^alfa2 on the uniform mesh
end                       
it       = zeros(1,2);
%
%  gamJ(:,(i-1)*s:i*s) = Jacobi-Fourier coefficients at step i, graded mesh
%  gamJc = Jacobi-Fourier coefficients at step i, uniform mesh
%
if flagalfa
   %%%% gamJ  --> gamJ1, gamJ2 
   gamJ1    = zeros(m1,nu*s);         % on the graded mesh                    
   gamJ2    = zeros(m2,nu*s);          
   %%%% gamJc --> gamJc1, gamJc2 
   gamJc1   = zeros(m1,(n-n1)*s);     % on the uniform mesh                   
   gamJc2   = zeros(m2,(n-n1)*s);
else
   gamJ1    = zeros(m,nu*s);          % on the graded mesh  
   gamJc1   = zeros(m,(n-n1)*s);      % on the uniform mesh
   gamJ2    = [];
   gamJc2   = [];
end   
indell   = 1:ell;
if flagalfa
   Im    = eye(m*s);
   ind1  = 1:m1;
   ind2  = m1+1:m;
else
   Im    = eye(m); 
   ind1  = 1:m;
   ind2  = [];
end   
%-------------------------------------------------------------------------%
flag_NaN = 0;
small    = 0.6;  % bound for using the blended or the fixed-point iteration
blend    = 0;    % number of blended steps
inds     = 0;
indk     = 0;
s1       = 1:s;
uno      = double(nu>n1); % 1 with also graded mesh, 0 if only uniform mesh
i        = 0;    % step counter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% graded mesh
for ii = 1:nu*uno
    i  = i+1;
    J0              = feval( fun, t(i), y(i,:), 1 ); 
    L               = min( norm(J0,1), norm(J0,inf) );
    if flagalfa
       blendflag    = L*fattb*max(halfa1(i),halfa2(i)) > small;
    else
       blendflag    = L*fattb*halfa1(i) > small;
    end   
    if blendflag
       if flagalfa 
          teta      = Im - [kron(halfa1(i)*X11,J0(ind1,ind1))  ...
                            kron(halfa2(i)*X12,J0(ind1,ind2)); ...
                            kron(halfa1(i)*X21,J0(ind2,ind1))  ...
                            kron(halfa2(i)*X22,J0(ind2,ind2))];
       else
          teta      = Im -halfa1(i)*gam*J0.'; 
       end   
       teta         = inv(teta);
       blend        = blend + 1;
    else
       teta         = [];
    end   
    tcor            = t(i) + hh(i)*x;
    tcor1           = [tcor; t(i+1)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% memory term and integration step
    if flagalfa
       fi0_1 = evalfi( y0(indell,ind1), tcor1, x1, ...
                                          Jj1(:,1:indk), gamJ1(:,1:inds) );   
       fi0_2 = evalfi( y0(indell,ind2), tcor1, x1, ...
                                          Jj2(:,1:indk), gamJ2(:,1:inds) );   
       fi0   = [fi0_1 fi0_2];
       %%%%%%%%%%%%%%%%%%%%%%%%
       [yi,gami,iti,flag_NaN] = step( fun, tcor, fi0, halfa1(i)*Is1, ...
                            halfa2(i)*Is2, PO1, PO2, fatt1_1*halfa1(i), ...
                            fatt1_2*halfa2(i), teta, L, m1 );
    else
       fi0   = evalfi( y0, tcor1, x1, Jj1(:,1:indk), gamJ1(:,1:inds) );
      %%%%%%%%%%%%%%%%%%%%%%%%
      [yi,gami,iti,flag_NaN] = step( fun, tcor, fi0, halfa1(i)*Is1, ...
                            [], PO1, [], fatt1_1*halfa1(i), ...
                            [], teta, L, m1, CC );
    end  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if blendflag 
       it(1)        = it(1) + iti;
    else   
       it(2)        = it(2) + iti;
    end   
    if flag_NaN
       fprintf( 'fhbvmdata: NaN at step %i, exit \n', i )            
       t = t(1:i);    y = y(1:i,:);                                      
       break
    end    
    y(i+1,:)         = yi;
    gamJ1(:,inds+s1) = ( fatt_1*halfa1(i) )*gami(:,ind1).';
    if flagalfa
       gamJ2(:,inds+s1) = ( fatt_2*halfa2(i) )*gami(:,ind2).';
    end   
    inds             = inds + s;
    indk             = indk + k+1;
end   
if flag_NaN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NaN: exit
   fprintf( 'NaN occurs at %g \n', t(i) ), return 
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data1.alfa  = data.alfa1;     %
data1.abd   = data.abd1;      %
data1.x     = data.x;         %
data1.x1    = data.x1;        %
data1.beta  = data.beta1;     %
data1.r     = data.r;         %
data1.nu    = data.nu;        %
data1.n1    = data.n1;        %
data1.Ps    = data.Ps1;       %
data1.glc   = data.glc;       %
data1.glb   = data.glb;       %
data1.s     = data.s;         %
data1.k     = data.k;         %
if flagalfa                   %
   data1.m     = data.m1;     %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
   data2.alfa  = data.alfa2;  %
   data2.abd   = data.abd2;   %
   data2.x     = data.x;      %
   data2.x1    = data.x1;     %
   data2.beta  = data.beta2;  %
   data2.r     = data.r;      %
   data2.nu    = data.nu;     %
   data2.n1    = data.n1;     %
   data2.Ps    = data.Ps2;    %
   data2.glc   = data.glc;    %
   data2.glb   = data.glb;    %
   data2.s     = data.s;      %
   data2.k     = data.k;      %
   data2.m     = data.m2;     %
else                          %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   data1.m  = m;              %
end                           %   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inds   = 0;
indk   = 0;
if flagalfa
   half12 = max(half1,half2);
else
   half12 = half1;
end   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uniform mesh
for j = n1+uno:n  
    i = i+1;       
    J0              = feval( fun, t(i), y(i,:), 1 ); 
    L               = min( norm(J0,1), norm(J0,inf) );
    blendflag       = L*fattb*half12 > small;
    if blendflag
       if flagalfa 
          teta      = Im - [kron(half1*X11,J0(ind1,ind1))  ...
                            kron(half2*X12,J0(ind1,ind2)); ...
                            kron(half1*X21,J0(ind2,ind1))  ...
                            kron(half2*X22,J0(ind2,ind2))];
       else
          teta      = Im - half1*gam*J0.';   
       end   
       teta         = inv(teta);
       blend        = blend + 1;
    else
       teta         = [];
    end   
    tcor            = t(i) + h*x;
    tcor1           = [tcor; t(i+1)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% memory term and integration step
    if flagalfa
      fi0_1 = evalfi( y0(indell,ind1), tcor1, x1, ...
                                        Jjc1(:,1:indk), gamJc1(:,1:inds) );
      fi0_2 = evalfi( y0(indell,ind2), tcor1, x1, ...
                                        Jjc2(:,1:indk), gamJc2(:,1:inds) );
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % sum of the contribution of the graded mesh to the memory term
      fi0_1 = fi0_1 + evalJjinu( j, gamJ1, data1 );
      fi0_2 = fi0_2 + evalJjinu( j, gamJ2, data2 );
      fi0   = [fi0_1 fi0_2];
      %%%%%%%%%%%%%%%%%%%%%%%%
      [yi,gami,iti,flag_NaN] = step( fun, tcor, fi0, half1*Is1, ...
          half2*Is2, PO1, PO2, fatt1_1*half1, fatt1_2*half2, teta, L, m1 );
    else
      fi0   = evalfi( y0, tcor1, x1, Jjc1(:,1:indk), gamJc1(:,1:inds) );
      fi0   = fi0 + evalJjinu( j, gamJ1, data1 );
      %%%%%%%%%%%%%%%%%%%%%%%%
      [yi,gami,iti,flag_NaN] = step( fun, tcor, fi0, half1*Is1, ...
                         [], PO1, [], fatt1_1*half1, [], teta, L, m1, CC );
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if blendflag 
       it(1)        = it(1) + iti;
    else   
       it(2)        = it(2) + iti;
    end   
    if flag_NaN
       fprintf( 'fhbvmdata: NaN at step %i, exit \n', i )           
       t = t(1:i);    y = y(1:i,:);                                  
       break
    end    
    y(i+1,:)          = yi;
    gamJc1(:,inds+s1) = ( fatt_1*half1 )*gami(:,ind1).';
    if flagalfa
       gamJc2(:,inds+s1) = ( fatt_2*half2 )*gami(:,ind2).';
    end   
    inds              = inds + s;
    indk              = indk + k+1;
end  
if flag_NaN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NaN: exit
   fprintf( 'NaN occurs at %g \n', t(i) )
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
   warning( 'findgam: blended iteration not A-convergent' )
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
for j   = 1:nu
    fi0 = fi0 + ( (t1j/fatj) )*y0(j,:);
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
    fi1  = fi1  + ( gamJ(:,indg+s1)*Jj(:,indJ+kk1) )';
    indg = indg - s;
    indJ = indJ + k1;
end
fi0 = fi0 + fi1;
return

%%%%%%%%%%%%%%%%%%%%%%%% The function below is specific for the hybrid mesh
function fi1 = evalJjinu( jj, gamJ, data ) 
%
%  Computes the contribute of the graded mesh to the memory term at step jj
%
%  gamJ are the Jacobi-Fourier coefficients relative to the initial 
%       graded mesh points.
%
%  data is the structure with the code data.
%
%                                                          Rel. 2025-09-24.
%%%%%%%%%%%%%%%%%%%%%
alfa  = data.alfa;  %
abd   = data.abd;   %
x     = data.x;     %
x1    = data.x1;    %
beta  = data.beta;  %
r     = data.r;     %
nu    = data.nu;    %
n1    = data.n1;    %
Ps    = data.Ps;    %
glc   = data.glc;   %
glb   = data.glb;   %
s     = data.s;     %
k     = data.k;     %
m     = data.m;     %
%%%%%%%%%%%%%%%%%%%%%
s1    = 1:s;
k1    = k+1;
% kk1   = 1:k1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  auxiliary integrals at step jj %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  and contribution to the memory term %%%%%%%
Jij   = zeros(s,k1);        
rnu1  = (r^nu-1)/n1;        
r1    = r-1;
ri1   = 1;    
fi1   = zeros(k1,m);
indg  = 0; 
for j = 1:nu
    den = ri1*r1;
    ind = 0;
    for i = 1:k1
        ind        = ind+1;
        M          = ( ( jj-1+x1(i) )*rnu1 + 1-ri1 )/den;
        Jij(:,ind) = Jjalfa1( M, abd, x, beta, alfa, Ps, glc, glb );       
    end   
    fi1  = fi1 + ( gamJ(:,indg+s1)*Jij )';
    indg = indg + s;
    ri1  = ri1*r;
end
return
%%%%%%%%%%%%%%%%%%%%%%%% The above function is specific for the hybrid mesh


function [y,gam,it,fNaN] = step( fun, t, fi0, Is1, Is2, PO1, PO2, ...
                                        ghalfa1, ghalfa2, teta, L, m1, CC )
%
%   Integration step for the local problem.
%                                                          Rel. 2025-10-03.
%
flagalfa = nargin==12; % 2 different alfas if true %%%%%%%
[k1,m]   = size(fi0);
k        = k1-1;
s        = size(PO1,1);
Y0       = fi0(1:k,:);
tol0     = 2*eps*(1+L);
tol      = 100*tol0;
itmax    = 1000*(m+1)*k1;
y        = Y0;
scal     = 1./( 1 + abs(Y0) );
err      = Inf;
gam      = zeros(s,m);   % Fourier coefficients
noblen   = isempty(teta); 
if flagalfa
   ind1  = 1:m1;
   ind2  = m1+1:m;
end
for it   = 1:itmax
    yold = y;
    f    = feval(fun,t,y);
    if noblen      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% fixed-point iteration              
       if flagalfa 
          gam  = [PO1*f(:,ind1) PO2*f(:,ind2)];     
       else
          gam  = PO1*f;
       end   
    else
       if flagalfa %%%%%%%%%%%%%%%%%%%%%%%%%%%% simplified Newton iteration
          eta  = gam - [PO1*f(:,ind1) PO2*f(:,ind2)]; 
          gam  = newtonit( gam, eta, teta, m1 );
       else        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% blended iteration
          eta  = gam - PO1*f;
          eta1 = CC*eta;
          gam  = gam -( eta1 +(eta-eta1)*teta )*teta;
       end    
    end        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% update the iterate
    if flagalfa    
       y = Y0 + [Is1*gam(:,ind1) Is2*gam(:,ind2)];
    else
       y = Y0 + Is1*gam;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check convergence
    err0 = err;
    err  = norm( (y-yold).*scal, inf );   
    if isnan(err) 
        disp( 'step: NaN' ), break
    elseif err<=tol0 
        break 
    elseif err>err0 && err0<=tol
        break
    end  
end                             
fNaN = any(any(isnan(y)));  %%%%%%%%%% a scalar flag is needed to check NaN
if ~fNaN 
   if err>tol*2 
      fprintf( 'step: step non converging: error = %g \n', err ) 
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% new approximation
   if flagalfa
      y = fi0(k1,:) + [ghalfa1*gam(1,ind1) ghalfa2*gam(1,ind2)]; 
   else
      y = fi0(k1,:) + ghalfa1*gam(1,:);
   end   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end   
return

function gam = newtonit( gam0, eta, teta, m1 )
%
%  Simplified Newton iteration.
%
%  Input:
%     gam0 : current approximation;
%     eta  : rhs of the Newton iteration;
%     teta : inverse of the simplified Jacobian;
%     m1   : equations with the first alpha.
%  Output:
%     gam  : new approximation.
%
%                                                          Rel. 2025-10-05.
[s,m] = size(gam0);
ind1  = 1:m1;
ind2  = m1+1:m;
m2    = m-m1;
eta1  = reshape( eta(:,ind1).',  s*m1, 1 );
eta2  = reshape( eta(:,ind2).',  s*m2, 1 );
gam01 = reshape( gam0(:,ind1).', s*m1, 1 );
gam02 = reshape( gam0(:,ind2).', s*m2, 1 );
gam   = [gam01; gam02] - teta*[eta1; eta2];
gam1  = reshape( gam(1:s*m1), m1, s ).';
gam2  = reshape( gam(s*m1+1:end), m2, s ).';
gam   = [gam1 gam2];                                     
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                         AUXILIARY FUNCTIONS                             %   
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,beta,PO,Is,Ps,Xs,abd] = jacobidat1( k, s, alfa, x1, beta1 )
%
%     [x,beta,PO,Is,Ps,Xs,abd] = jacobidat1(k,s,alfa)
%   
%                                   for Gauss-Jacobi abscissae and weights,
%   or
%
%     [x,beta,PO,Is,Ps,Xs,abd] = jacobidat1(k,s,alfa,x1,beta1)
%
%                                 for Jacobi-Pineiro abscissae and weights, 
%
%  along with auxiliary matrices for FHBVM(k,s).
%
%  INPUT:
%       - k,s  : parameters of the FHBVM(k,s) method, 1<=s<=k;
%       - alfa : index of the fractional derivative (alfa>0 and not in IN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       - x1    : abscissae where to compute Ps and Is;
%       - beta1 : wights quadrature, with unit sum;
%       if they are not given in input, they are computed as the
%       Gauss-Jacobi (alfa-1,0) ones.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  OUTPUT:
%       - x    : k Gauss-Jacobi/Jacobi-Pineiro abscissae in [0,1], relative
%                to the weight function    w(x) = alfa*(1-x)^alfa;
%       - beta : corresponding quadrature weights; ( sum(beta)==1 )
%       - PO   : PO = P_s^T * Omega, with P_s(i,j) = P_{j-1}(x_i)
%                                    and Omega = diag(beta);
%       - Is   : matrix I_s(i,j) = I^alfa P_{j-1}(x_i); % i = 1,...,k,
%       - Ps   : matrix P_s(i,j) = P_{j-1}(x_i);        % j = 1,...,s.
%       - Xs   : = PO*Is;
%       - abd  : coefficients [a b d] for the recurrence:
%               P_{j+1}(c)=(a(j)*c-b(j))*P_j(c)-d(j)*P_{j-1}(c), j=0,..s-2,
%               P_0(c)==1, P_{-1}(c)==0.
%
%                                                          Rel. 2025-10-05.
%
if (s>k) || (s<1) || (alfa<=0) || (alfa==round(alfa)) 
    error( 'jacobidat1: wrong parameters' ) 
end
ab   = r_jacobi01(k,alfa-1,0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==3
   xw   = gauss(k,ab);
   x    = xw(:,1);
   beta = xw(:,2);
   beta = beta*alfa; % so that sum(beta)==1;
   flag = 1;
   if nargout<=2, return, end
else
   x    = x1(:);
   beta = beta1(:);  % sum(beta1)==1;
   flag = 0;
end   
if abs(sum(beta)-1)>10*k*eps % check that the Jacobi weights are accurate %
    disp( 'jacobidat1: are weights accurate? please check' ) 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ps, PO
Ps      = zeros(k,k+1);
Ps(:,1) = 1;
Ps(:,2) = (x-ab(1,1)).*Ps(:,1);
for j  = 2:k
    Ps(:,j+1) = (x-ab(j,1)).*Ps(:,j) - ab(j,2)*Ps(:,j-1);
end
if flag
   err = norm(Ps(:,k+1),inf);  
   if err>k*eps   % check that the Jacobi abscissae are accurate %
      disp( 'jacobidat1: are abscissae accurate? please check' )   
   end
end   
Ps    = Ps(:,1:s);
% polynomials are scaled in order to be orthonormal
scal  = 1./sqrt(diag(Ps.'*diag(beta)*Ps));
Ps    = Ps*diag(scal);
PO    = Ps.'*diag(beta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Is, Xs
Is    = zeros(k,s);
for i = 1:k    %, disp([i k])
    Is(i,:) = makeint( ab(:,1), ab(:,2), alfa, x(i), s, x, beta );  
end    
Is    = Is*diag(scal);
Xs    = PO*Is;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% abd, validated 2023.08.21
abd   = zeros(s-1,3);
dx    = diag(x);
for j = 1:s-1
    aj = 1/( PO(j+1,:)*diag(x)*Ps(:,j) );
    bj = aj*PO(j,:)*dx*Ps(:,j);
    if j==1
       dj    = 0;
    else
       dj    = aj/abd(j-1,1);
    end    
    abd(j,:) = [aj bj dj];
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function [Ps1,Ps2,x,beta] = makePs1( abd1, abd2, k )
%
%   Compute the k Gauss-Legendre nodes (x) and weights (beta), along with
%   the Jacobi matrices evaluated at the Gauss-abscissae (Ps1 and Ps2).
%   abd1 and abd2 contain the coefficients of the Jacobi recursion.
%
%                                                          Rel. 2025-10-02.
%   REMARK: order of the Gauss-Legendre quadrature = 2*k
%
a1       = abd1(:,1);
b1       = abd1(:,2);
d1       = abd1(:,3);
flag     = ~isempty(abd2);  
if flag
   a2    = abd2(:,1);
   b2    = abd2(:,2);
   d2    = abd2(:,3);
else
   Ps2   = []; 
end   
s        = size(abd1,1)+1;
[x,beta] = GL(k);  % Gauss-Legendre abscissae and weights
x        = x(:);
beta     = beta(:);
Ps1      = zeros(k,s);
Ps1(:,1) = 1;
if flag
   Ps2      = zeros(k,s);
   Ps2(:,1) = 1;
end   
if s>1 
   Ps1(:,2)     = a1(1)*x-b1(1);
   if flag
      Ps2(:,2)  = a2(1)*x-b2(1);
   end   
   for j = 2:s-1
       Ps1(:,j+1)    = (a1(j)*x-b1(j)).*Ps1(:,j) -d1(j)*Ps1(:,j-1);
       if flag
          Ps2(:,j+1) = (a2(j)*x-b2(j)).*Ps2(:,j) -d2(j)*Ps2(:,j-1);
       end   
   end
end   
return

function Ialfa = Jjalfa1( M, abd, cJ, wJ, alfa, Ps, x, beta )
%
%   Ialfa = Jjalfa1( M, abd, cJ, wJ, alfa, Ps, x, beta )  
%
% Computes the integrals:
% 
%       J_j^alfa(M) = int_0^1(M-c)^(alfa-1)P_j(c)dc, j = 0...s-1,
%
% with  {P_j(c)} satisfying the recurrence:
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
%
%                                                          Rel. 2025-02-07.
%
Mmax = 1.1;  
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
% {P_j(c)} satisfying the recurrence:
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
wit     = wi(:).';
Ps(:,1) = 1;
if s>1
   xx      = x*ci;
   Ps(:,2) = a(1)*xx - b(1);
end   
for j = 2:s-1
    Ps(:,j+1) = ( a(j)*xx - b(j) ).*Ps(:,j) -d(j)*Ps(:,j-1);
end
Ialfa   = x^alfa*( wit*Ps );
Ps(:,1) = 1;
if s>1
   xx      = 1 + (x-1)*ci;
   Ps(:,2) = a(1)*xx - b(1);
end   
for j = 2:s-1
    Ps(:,j+1) = ( a(j)*xx - b(j) ).*Ps(:,j) -d(j)*Ps(:,j-1);
end
Ialfa = Ialfa - (x-1)^alfa*( wit*Ps );
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
%
%   ab = r_jacobi01(N,a,b)
%
%   Recurrence coefficients for monic shifted Jacobi polynomials.
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
%
%                                                       Restyle 2025-09-29.
%
if nargin<2, a = 0; end;  if nargin<3, b = a; end
if((N<=0)||(a<=-1)||(b<=-1))
    error( 'r_jacobi01: parameter(s) out of range' ) 
end
cd = r_jacobi(N,a,b);
n  = 1:N;  ab(n,1) = (1+cd(n,1))./2;  
ab(1,2) = cd(1,2)/2^(a+b+1);
n = 2:N;  ab(n,2) = cd(n,2)./4;
return

function xw = gauss(N,ab)
%
%   xw = gauss(N,ab)                    Gauss quadrature rule.
%
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
%
%                                                       Restyle 2025-09-29.
%
N0 = size(ab,1); 
if N0<N, error( 'gauss: input array ab too short' ), end
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
%
%   ab = r_jacobi(N,a,b)
%
%   Recurrence coefficients for monic Jacobi polynomials.
%   
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
%
%                                                       Restyle 2025-09-29.
%
if nargin<2, a = 0; end  
if nargin<3, b = a; end
if((N<=0)||(a<=-1)||(b<=-1)) 
    error( 'r_jacobi: parameter(s) out of range' ) 
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
% [x,w,y] = lgwt01(N)
%
% This script computes the Legendre-Gauss nodes (x) and weights 
% (w) on the interval [0,1] with truncation order N.
% Also the nodes (y) on [-1,1] are computed.
% Adapted from the file written by Greg von Winckel-25/02/2004.
% 
%                                                          Rel. 2024-07-17.
%
a  = 0; b = 1;  %%%%%%%%%%%%%%%%%%%%% normalization on the interval [0,1] %
N  = N-1;
N1 = N+1; N2 = N+2;
xu = linspace(-1,1,N1)';
% Initial guess %
y  = cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
L  = zeros(N1,N2);   % Gauss-Legendre Vandermonde Matrix %
% Compute the zeros of the N+1 Legendre Polynomial using %
% the recursion relation and the Newton-Raphson method   %
y0 = 2;                    % Iterate until new points are uniformly %
while max(abs(y-y0))>eps   %  within epsilon of old points          %    
    L(:,1) = 1;  L(:,2) = y; 
    for k = 2:N1
        L(:,k+1) = ( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
    Lp = (N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    y0 = y;  y = y0-L(:,N2)./Lp;
end
x = (a*(1-y)+b*(1+y))/2;      % tranformation from [-1,1] to [a,b]==[0,1] %    
w = (b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;             % Compute the weights %
w = ( w + w(end:-1:1) )/2;                          % and symmetrize them %
return
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

function [x,w1,w2,ier] = myJacobiPineiro( k, alpha1, alpha2 )
%
%   [x,w1,w2,ier] = myJacobiPineiro( k, alpha1, alpha2 )
%
%  Computes the Jacobi-Pineiro abscissae and the weights of the quadratures
%
%   \sum_{i=1}^k w1_i f(x_i)    
%   \sum_{i=1}^k w2_i f(x_i)
%
%   which are, respectively, exact for the integrals
%
%   alpha1\int_0^1 (1-x)^{alpha1-1} f(x) dx
%   alpha2\int_0^1 (1-x)^{alpha2-1} f(x) dx
%
%   with alpha1 \ne alpha2, alpha1,alpha2 \in(0,1),
%
%   for all polynomials f(x) of degree no larger than [1.5*k]-1,
%   with [] the integer part of the argument.
%
%   Input:
%     k - better if even, so that [1.5*k] = 1.5*k = 3*(k/2);
%     alpha1, alpha2 - they need to be positive and their difference has 
%                      not to be an integer. 
%
%   Output:
%     x - abscissae;
%     w1, w2 - weights of the two formulae;
%     ier - error code, if different from 0.
%
%   Adapted from:
%   T.Laudadio, N.Mastronardi, W.Van Assche, P.Van Dooren
%   A Matlab package computing simultaneous Gaussian quadrature rules for
%   multiple orthogonal polynomials
%   J. Comput. Appl. math. 451 (2024) 116109
%   https://doi.org/10.1016/j.cam.2024.116109
%
%                                                          Rel. 2025-10-05.
%
if nargin<3
    error( 'number of input parameters wrong' )
elseif abs(alpha1-alpha2)==fix(abs(alpha1-alpha2)) || min(alpha1,alpha2)<=0
    error( 'wrong alphas' )
end
[b,c,d,F] = MOP1(k,[0 alpha1-1 alpha2-1]);
[x,w1,w2,ier] = GaussMOP1( b, c, d, k, F );
w1 = w1*alpha1; % => sum(w1)==1;
w2 = w2*alpha2; % => sum(w2)==1;
x  = 1-x;       % x -> 1-x;
return 

%%%%%%%%%% auxiliary functions adapted from the above reference %%%%%%%%%%%

function [b,c,d,F] = MOP1(n,v)
%
% [b,c,d,F] = MOP1(n,v)
%
% computations of the recurrence relation coefficients
% associated with the multiple Jacobi-Pineiro orthogonal polynomial.
%
%                                                       Restyle 2025-09-29.
%
a0 = v(1);
a1 = v(2);
a2 = v(3);
[b,c,d,F] = Multiple_Jacobi_Pineiro(a0,a1,a2,n);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[bb,cc,dd,F] = Multiple_Jacobi_Pineiro(a0,a1,a2,n)
%
% [bb,cc,dd,F] = Multiple_Jacobi_Pineiro(a0,a1,a2,n)
%
% computation of the coefficients of the recurrence relation associated 
% with the multiple Jacobi-Pineiro polynomials of first kind
%  
% Reference: 
%   W. Van Assche, E. Coussement 
%   Some classical multiple orthogonal polynomials 
%   J. Comput. Appl. Math. 127 (2001) 317-347.
%
%                                                       Restyle 2025-09-29.
%
m  = floor(n/2);
%%%%%%%%%%%%%%%%%%%%
bb = zeros(1,2*m+2);
cc = bb;
dd = bb; 
%%%%%%%%%%%%%%%%%%%%
for i = 0:m    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bb(2*i+1) = (36*i^4+(48*a0+28*a1+20*a2+38)*i^3+ ...
      (21*a0^2+8*a1^2+4*a2^2+30*a0*a1+18*a0*a2+15*a1*a2+39*a0+ ...
      19*a1+19*a2+9)*i^2 ...
      +(3*a0^3+10*a0^2*a1+4*a0^2*a2+6*a0*a1^2+2*a0*a2^2+11*a0*a1*a2+ ...
      5*a1^2*a2+3*a1*a2^2+12*a0^2+3*a1^2+3*a2^2+13*a0*a1+ ...
      13*a0*a2+8*a1*a2+6*a0+3*a1+3*a2)*i+ ...
      a0^2+a0*a1+a2*a1^2+2*a2*a1^2*a0+2*a0^2*a1+a1^2*a0+a2^2*a0+ ...
      a2^2*a1+a0^3*a1+a0^2*a1^2+a2^2*a0*a1+a2^2*a1^2+2*a2*a0^2*a1 ...
      +3*a2*a1*a0+2*a2*a0^2+a1*a2+a0^3+a0*a2)/ ...
      ((3*i+a0+a2)*(3*i+a0+a1)*(3*i+a0+a2+1)*(3*i+a0+a1+2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bb(2*i+2) = (36*i^4+(48*a0+20*a1+28*a2+106)*i^3+ ...
      (21*a0^2+4*a1^2+8*a2^2+18*a0*a1+30*a0*a2+15*a1*a2+ ...
      105*a0+41*a1+65*a2+111)*i^2 ...
      +(3*a0^3+4*a0^2*a1+10*a0^2*a2+2*a0*a1^2+6*a0*a2^2+11*a0*a1*a2+ ...
      3*a1^2*a2+5*a1*a2^2 ...
      +30*a0^2+5*a1^2+13*a2^2+23*a0*a1+47*a0*a2+22*a1*a2+72*a0+25*a1+ ...
      49*a2+48)*i ...
      +18*a0*a2+8*a2*a0^2+4*a1+4*a2^2*a1+8*a1*a2+2*a0^3+5*a2^2*a0+ ...
      8*a2*a1*a0+12*a2 ...
      +7+15*a0+a2^2*a1^2+10*a0^2+6*a0*a1+2*a2*a1^2+2*a0^2*a1+a1^2*a0+ ...
      5*a2^2+a2*a0^3 ...
      +a2^2*a0^2+a1^2+a2*a1^2*a0+2*a2*a0^2*a1+2*a2^2*a0*a1)/ ...
      ((3*i+a0+a2+1)*(3*i+a0+a1+2)*(3*i+a0+a2+3)*(3*i+a0+a1+3));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cc(2*i+1) = (i*(2*i+a0)*(2*i+a0+a1)*(2*i+a0+a2)*(54*i^4+ ...
      (63*a0+45*a1+45*a2)*i^3 ...
      +(24*a0^2+8*a1^2+8*a2^2+42*a0*a1+42*a0*a2+44*a1*a2-8)*i^2+ ...
      (3*a0^3+a1^3+a2^3+12*a0^2*a1+12*a0^2*a2+3*a0*a1^2+ ...
      3*a0*a2^2+33*a0*a1*a2 ...
      +8*a1^2*a2+8*a1*a2^2-3*a0-4*a1-4*a2)*i+a0^3*a1+a0^3*a2+...
      6*a0^2*a1*a2+a1^3*a2+a1*a2^3+3*a0*a1^2*a2+3*a0*a1*a2^2-a0*a1- ...
      a0*a2-2*a1*a2 ...
      ))/((3*i+a0+a1+1)*(3*i+a0+a2+1)*(3*i+a0+a1)^2*(3*i+a0+a2)^2* ...
      (3*i+a0+a1-1)*(3*i+a0+a2-1)); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cc(2*i+2) = ((2*i+a0+1)*(2*i+a0+a1+1)*(2*i+a0+a2+1)*(54*i^5+ ...
      (63*a0+45*a1+45*a2+135)*i^4 ...
      +(24*a0^2+8*a1^2+8*a2^2+42*a0*a1+42*a0*a2+44*a1*a2+126*a0+76*a1+ ...
      104*a2+120)*i^3+(3*a0^3+a1^3+a2^3+12*a0^2*a1+12*a0^2*a2+ ...
      3*a0*a1^2+3*a0*a2^2+33*a0*a1*a2 ...
      +8*a1^2*a2+8*a1*a2^2+36*a0^2+5*a1^2+19*a2^2+54*a0*a1+72*a0*a2+ ...
      66*a1*a2+87*a0+39*a1+81*a2+45)*i^2 ...
      +(a0^3*a1+a0^3*a2+6*a0^2*a1*a2+a1^3*a2+a1*a2^3+3*a0*a1^2*a2+3* ...
      a0*a1*a2^2+3*a0^3+2*a2^3+ ...
      12*a0^2*a1+12*a0^2*a2+6*a0*a2^2+33*a0*a1*a2+5*a1^2*a2+11*a1*a2^2+ ...
      18*a0^2+20*a0*a1+38*a0*a2+14*a2^2+26*a1*a2+24*a0+6*a1+24*a2+6)*i+ ...
      a0^3*a1+3*a0^2*a1*a2+3*a0*a1*a2^2+a1*a2^3+a0^3+a2^3+3*a0^2*a1+ ...
      3*a0^2*a2+6*a0*a1*a2+3*a0*a2^2+3*a1*a2^2+3*a0^2+3*a2^2+2*a0*a1+ ...
      6*a0*a2+2*a1*a2+2*a0+2*a2)) ...
      /((3*i+a0+a1+3)*(3*i+a0+a2+2)*(3*i+a0+a1+2)^2*(3*i+a0+a2+1)^2* ...
      (3*i+a0+a1+1)*(3*i+a0+a2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dd(2*i+1) = (i*(2*i+a0)*(2*i+a0-1)*(2*i+a0+a1)*(2*i+a0+a1-1)* ...
      (2*i+a0+a2)*(2*i+a0+a2-1)*(i+a1)*(i+a1-a2))/ ...
      ((3*i+a0+a1+1)*(3*i+a0+a1)^2*(3*i+a0+a2)*(3*i-1+a0+a1)^2* ...
      (3*i+a0+a2-1)*(3*i+a0+a1-2)*(3*i+a0+a2-2));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    dd(2*i+2) = (i*(2*i+a0+1)*(2*i+a0)*(2*i+a0+a1)*(2*i+a0+a1+1)* ...
      (2*i+a0+a2+1)*(2*i+a0+a2)*(i+a2)*(i+a2-a1))/ ...
      ((3*i+a0+a1+2)*(3*i+a0+a2+2)*(3*i+a0+a1+1)*(3*i+1+a0+a2)^2* ...
      (3*i+a0+a1)*(3*i+a0+a2)^2*(3*i+a0+a2-1));
end
dd(1)  = 0;
dd(2)  = 0;
cc(1)  = 0;
bb(1)  = (1+a1)/(2+a0+a1);
cc(2)  = (1+a0)*(1+a1)/((3+a0+a1)*(2+a0+a1)^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F(1,1) = gamma(1+a0)*gamma(1+a1)/gamma(2+a0+a1);
F(2,1) = gamma(1+a0)*gamma(1+a2)/gamma(2+a0+a2);
F(2,2) = ((1+a2)-(2+a0+a2)*bb(1))*gamma(1+a0)*gamma(1+a2)/gamma(3+a0+a2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[x,w1,w2,ier] = GaussMOP1(b,c,d,n,F)
%
% [x,w1,w2,ier] = GaussMOP1(b,c,d,n,F)
%
% computation of the nodes oand weights of the simultaneous Gaussian
% quadrature rule associated to multiple orthogonal polynomials
%
% Input:
%   b : diag({H}_{n});
%   c : diag({H}_{n});
%   d : diag({H}_{n},-2);
%   n : size of \hat{H}_n
%   F : 2x2  matrix used for computing the weights
%
% Output:
%   x, w1, w2 : nodes and weights of the simultaneous Gauss rule
%   ier       : ier is set to zero for normal return, ier is set to j if 
%          the j-th eigenvalue has not been determined after 30 iterations.
%
%                                                       Restyle 2025-10-04.
%

% balancing the Hessenberg matrix
a = ones(1,n+1); [as,bs,cs,ds,S] = DScaleS2(a,b,c,d,n+1);
%
% transformation of the Hessenberg matrix into a similar nonsymmetric 
% tridiagonal matrix by means of elementary transformations
[~,as1,bs1,cs1] = tridEHbackwardV(as,bs,cs,ds,n);
%
% transformation of the nonsymmetric tridiagonal matrix into a similar 
% symmetric tridiagonal matrix
[bs2,cs2] = DScaleSV2(as1,bs1,cs1,n+1);
%
% computation of the eigenvalues of the similar symmetric tridiagonal 
% matrix
%[x,e1,z1,ierr]= gausq2(n, bs2(1:n), cs2(2:n));
[x,~,~,~]= gausq2(n, bs2(1:n), cs2(2:n));  x = sort(x);
%
% computation of the nodes and weights
% by means of a modification of the Ehrlich–Aberth method
tol = n^2*eps; [x,w1,w2,ier] = EA_method(as,bs,cs,ds,x,tol,F,S,n);
%
return 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a,b,c,d,S] = DScaleS2(a,b,c,d,n)
%
%   [a,b,c,d,S] = DScaleS2(a,b,c,d,n)
%
%   Diagonal scaling of the  banded Hessenberg matrix
%   has as input a lower Hessenberg matrix H with bandwidth 4 
%
%   Input:
%     a = diag(H,1);
%     b = diag(H);
%     c = diag(H,-1);
%     d = diag(H,-2);
%
%   and returns a scaled version Hs = S\H*S, 
%   where S is a diagonal matrix, of the input vectors.
%
%                                                       Restyle 2025-09-29.
%
S(1) = 1; 
S(2) = sqrt(c(2)); 
t1   = sqrt(a(1)/c(2));
a(1) = a(1)/t1;
c(2) = a(1);
for i = 2:n-1
    t2     = sqrt(a(i)/c(i+1));
    d(i+1) = d(i+1)*t1*t2;
    a(i)   = sqrt(a(i)*c(i+1));
    c(i+1) = a(i); 
    t1     = t2;
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GEL,a,b,c] = tridEHbackwardV(a,b,c,d,n)
%
%   [GEL,a,b,c] = tridEHbackwardV(a,b,c,d,n)
%
%   reduction of the lower Hessenberg matrix  \hat{H}_n
%   into the similar  similar tridiagonal one \hat{T}_n
%
%   Input: 
%
%     a = diag(\hat{H}_n,1);
%     b = diag(\hat{H}_n);
%     c = diag(\hat{H}_n,-1);
%     d = diag(\hat{H}_n,-2);
%     n : size of \hat{H}_n;
%
%   Output:
%
%     a = diag(\hat{T}_n,1);
%     b = diag(\hat{T}_n);
%     c = diag(\hat{T}_n,-1);
%
%                                                       Restyle 2025-09-29.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%
last = fix(n*(n-2)/4);
GEL  = zeros(last,3); 
%%%%%%%%%%%%%%%%%%%%%%%%%%
ind  = 0;
for i = n:-1:3
    %
    hh1    = d(i)/c(i);
    ind    = ind+1;
    %
    c(i-1) = c(i-1)-hh1*b(i-1);
    b(i-2) = b(i-2)-hh1*a(i-2);
    %
    d(i)   = 0;
    t      = d(i-2)*hh1;
    d(i-1) = d(i-1)+hh1*c(i-2);
    c(i-1) = c(i-1)+hh1*b(i-2);
    b(i-1) = b(i-1)+hh1*a(i-2);
    %%%%%%%%%%%%%%%%
    for j = i-1:-2:4
        %
        hh1    = t/d(j);
        d(j-1) = d(j-1)-hh1*c(j-1);
        c(j-2) = c(j-2)-hh1*b(j-2);
        b(j-3) = b(j-3)-hh1*a(j-3);
        %
        ind    = ind+1;
        GEL(ind,:) = [hh1,j,j-3];
        %
        t      = d(j-3)*hh1;
        d(j-2) = d(j-2)+hh1*c(j-3);
        c(j-2) = c(j-2)+hh1*b(j-3);
        b(j-2) = b(j-2)+hh1*a(j-3);
        %
    end
    %%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%
GEL(ind+1:last,:) = [];
%%%%%%%%%%%%%%%%%%%%%%%
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [b,c] = DScaleSV2(a,b,c,n)
%
%   [b,c] = DScaleSV2(a,b,c,n)
%
%   transformation of the unsymmetric tridiagonal matrix T
%   into  a similar  symmetric tridiagonal one
%
%   Input:
%
%     a = diag(T,1);
%     b = diag(T);
%     c = diag(T,-1);
%     n : size of the matrix.
%
%   Output:
%
%     a = diag(T,1);
%     b = diag(T);
%     c = diag(T,-1).   
%
%                                                       Restyle 2025-09-29.
%
for i = 1:n-1
    a(i)   = sqrt(a(i)*c(i+1));
    c(i+1) = a(i);
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w1,w2,ier] = EA_method(a,b,c,d,x,tol,D,S,n)
%
%  [x,w1,w2,ier] = EA_method(a,b,c,d,x,tol,D,S,n)
%
%  Ehrlich-Aberth Method for computing the zeros and the left and right 
%  eigenvector of the Hessenberg matrix associated to multiple orthogonal 
%  polynomials and the nodes and weights of the simultaneous Gauss rule
%
%  Input:
%
%    a   = diag( \hat{H}_{n,n+1},1);
%    b   = diag( \hat{H}_{n,n+1});
%    c   = diag( \hat{H}_{n,n+1},-1);
%    d   = diag( \hat{H}_{n,n+1},-2);
%    x   : vector of the approximation of eigenvalues of \hat{H}_n;
%    tol : tolerance in computing the eigenvalues;
%    D   : 2x2 matrix used for computing the weights;
%    S   : first two diagonal entries of S_n;
%    n   : size of \hat{H}_n.
%
% Output:
%
%   x      : nodes of the simultaneous Gauss rule
%   w1, w2 : weights of the simultaneous Gauss rule
%   ier    : error flag
%            ier is set to zero for normal return,
%            ier is set to j if the j-th eigenvalue has not been
%            determined after 30 iterations.
%
%                                                       Restyle 2025-09-29.
%
m     = n;
n     = n+1;
convj = zeros(m,1);
iter  = 0;
ier   = 0;
iterj = zeros(m,1);
%%%%%%%%%%%%%%%%%%%
w1    = zeros(1,m);
w2    = zeros(1,m);
%%%%%%%%%%%%%%%%%%%
while sum(convj)<m && ier==0  %% ier is an error flag
    iter = iter+1;
    %%%%%%%%%%%
    for j = 1:m
        x0 = x(j);
        if convj(j) == 0
            % computation of the new approximation of the eigenvalue x(j)
            [p0,p1] = one_step_Newton_nullspaceV(a,b,c,d,x0,n);
            p01     = p0(n)/p1(m);
            %
            sum1  = 0;
            for k = 1:m
                if k~=j
                    sum1 = sum1+1/(x(j)-x(k));
                end
            end
            t     = x(j)-p01/(1-p01*sum1);
            x(j)  = t;
            if abs(p01)< tol || abs(p0(n)) < tol
               convj(j) = 1;   
              % computation of the right eigenvector v associated with x(j)
               v = p0(1:m);
              % computation of the left eigenvector u associated with x(j)
                [u,~,~,~,~] = compute_left_V(x(j),a,b,c,d,m);
              % computation of the scalar product u^T v
                uv    = dot(u,v);
              % computation of the weights of the simultaneous Gauss
              % quadrature rule
                w1(j) = v(1)*(D(1,1)*u(1))/uv;
                w2(j) = S(1)*v(1)*(D(2,1)*u(1)/S(1)+D(2,2)*u(2)/S(2))/uv;
              %
                iterj(j) = iter;
            end
        end
    end
    %%%%%%%%%%%%%%
    if iter==30  %
        ier = j; %  for the stop criterion 
    end          %
    %%%%%%%%%%%%%%
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [p,p10] = one_step_Newton_nullspaceV(a,b,c,d,lambda,n)
%
%  [p,p10] = one_step_Newton_nullspaceV(a,b,c,d,lambda,n)
%
%     computation of the sequence of  MOPs and their derivatives evaluated
%     in lambda
%
% Input:
%
%   a      = diag(\hat{H}_{n,n+1},1);
%   b      = diag(\hat{H}_{n,n+1});
%   c      = diag(\hat{H}_{n,n+1},-1);
%   d      = diag(\hat{H}_{n,n+1},-2)
%   lambda : approximation of an eigenvalue of  \hat{H}_n
%   n      : size of \hat{H}_n.
%
% Output:
%
%   p   : vector [p_0(lambda), p_1(lambda),\ldots,p_n(lambda)]^T;
%   p10 : vector [p_1(lambda), p_2(lambda),\ldots,p_n(lambda)]^T.
%
%                                                       Restyle 2025-09-29.
%
m  = n;
n  = n+1;
b  = b-lambda;
a0 = a;
b0 = b;
c0 = c;
d0 = d;
p  = zeros(n,1);
%
% transformation of the nx(n+1) Hessenberg matrix into a lower tridiagonal
% one by means of the multiplication of n Givens rotations to the right
% in order to annihilate the first superdiagonal of the Hessenberg matrix
%%%%%%%%%%%%%%%%%
CS1 = zeros(m,2);
%%%%%%%%%%%%%%%%%
for i = 1:m-3
    G        = givens_GV(b(i),a(i)); 
    CS1(i,1) = G(1,1); 
    CS1(i,2) = G(1,2);
    V1       = [b(i),a(i)]*G';
    b(i)     = V1(1); 
    a(i)     = V1(2);
    V1       = [c(i+1),b(i+1)]*G';
    c(i+1)   = V1(1);
    b(i+1)   = V1(2);
    V1       = [d(i+2),c(i+2)]*G';
    d(i+2)   = V1(1);
    c(i+2)   = V1(2); 
    V1       = [0,d(i+3)]*G';
    d(i+3)   = V1(2);
end
i = m-2; %%%%%%%%%%%%%%%%%%%%%%%%%
  G        = givens_GV(b(i),a(i));
  V1       = [b(i),a(i)]*G';
  b(i)     = V1(1);
  a(i)     = V1(2);
  V1       = [c(i+1),b(i+1)]*G';
  c(i+1)   = V1(1);b(i+1)=V1(2);
  V1       = [d(i+2),c(i+2)]*G';
  d(i+2)   = V1(1);
  c(i+2)   = V1(2);
  CS1(i,1) = G(1,1); 
  CS1(i,2) = G(1,2);
i = m-1; %%%%%%%%%%%%%%%%%%%%%%%%%
  G        = givens_GV(b(i),a(i));
  V1       = [b(i),a(i)]*G';
  b(i)     = V1(1);
  a(i)     = V1(2);
  V1       = [c(i+1),b(i+1)]*G';
  c(i+1)   = V1(1);b(i+1)=V1(2);
  CS1(i,1) = G(1,1); 
  CS1(i,2) = G(1,2);
i = m;   %%%%%%%%%%%%%%%%%%%%%%%%%
  G        = givens_GV(b(i),a(i));
  V1       = [b(i),a(i)]*G';
  b(i)     = V1(1);
  a(i)     = V1(2);
  CS1(i,1) = G(1,1); 
  CS1(i,2) = G(1,2);
%
% computation of the vector 
%              [p_0(lambda); p_1(lambda);...;p_{n-1}(lambda);p_{n}(lambda)]
%
p(n) = CS1(n-1,1);
t    = 1;
for i = n-1:-1:2
    t    = -t*CS1(i,2);
    p(i) = CS1(i-1,1)*t;
end
p(1) = -t*CS1(1,2);
%
% computation of the vector 
%                        [p'_1(lambda);...;p'_{n-1}(lambda);p'_{n}(lambda)]
%
p10    = zeros(m,1);
p10(1) = p(1)/a0(1);
p10(2) = (p(2)-b0(2)*p10(1))/a0(2);
p10(3) = (p(3)-c0(3)*p10(1)-b0(3)*p10(2))/a0(3);
for i  = 4:m
    p10(i) = (p(i)-d0(i)*p10(i-3)-c0(i)*p10(i-2)-b0(i)*p10(i-1))/a0(i);
end
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[u1,a,b,c,d] = compute_left_V(lambda,a,b,c,d,n)
%
%  [u1,a,b,c,d] = compute_left_V(lambda,a,b,c,d,n)
%
%  computation of the left eigenvector u associated with the eigenvalue 
%  lambda of the n x n Hessenberg matrix by means of the multiplication 
%  of  n-1 Givens rotations to the left in order to annihilate the first 
%  superdiagonal of the Hessenberg matrix
%
% Input:
%
%   lambda : approximation of an eigenvalue of  \hat{H}_n
%   a      = diag(\hat{H}_{n,n+1},1);
%   b      = diag(\hat{H}_{n,n+1});
%   c      = diag(\hat{H}_{n,n+1},-1);
%   d      = diag(\hat{H}_{n,n+1},-2);
%   n      : size of \hat{H}_n.
%
% Output:
%
%   u1     : left eigenvector associated to lambda.
%
%                                                       Restyle 2025-09-29.
%
b = b-lambda;
for i = n-1:-1:3
    G       = givens_GV(b(i+1),a(i));
    cs(i,1) = G(1,1); 
    cs(i,2) = G(1,2);
    V       = G'*[d(i) c(i) b(i) a(i); 0 d(i+1) c(i+1) b(i+1)];
    d(i)    = V(1,1); 
    c(i)    = V(1,2);
    b(i)    = V(1,3);
    a(i)    = V(1,4);
    d(i+1)  = V(2,2); 
    c(i+1)  = V(2,3); 
    b(i+1)  = V(2,4);
end
i = 2; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G       = givens_GV(b(i+1),a(i));
    cs(i,1) = G(1,1); 
    cs(i,2) = G(1,2);
    V       = G'*[c(i) b(i) a(i); d(i+1) c(i+1) b(i+1)];
    c(i)    = V(1,1);
    b(i)    = V(1,2);
    a(i)    = V(1,3);
    d(i+1)  = V(2,1); 
    c(i+1)  = V(2,2); 
    b(i+1)  = V(2,3);
i = 1;  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    G       = givens_GV(b(i+1),a(i));
    cs(i,1) = G(1,1); 
    cs(i,2) = G(1,2);
    V       = G'*[b(i) a(i);  c(i+1) b(i+1)];
    b(i)    = V(1,1);
    a(i)    = V(1,2);
    c(i+1)  = V(2,1); 
    b(i+1)  = V(2,2);
    t       = 1; 
    u1      = zeros(n,1);
    for i = 1:n-1
        u1(i) =  t*cs(i,1);
        t     = -t*cs(i,2);
    end
    u1(n) = t;
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d,e,z,ierr] = gausq2(n, d, e)
%
%     [d,e,z,ierr] = gausq2(n, d, e)
%
%     this subroutine is a translation of an algol procedure,
%     num. math. 12, 377-383(1968) by martin and wilkinson,
%     as modified in num. math. 15, 450(1970) by dubrulle.
%     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
%     this is a modified version of the 'eispack' routine imtql2.
%
%     this subroutine finds the eigenvalues and first components of the
%     eigenvectors of a symmetric tridiagonal matrix by the implicit ql
%     method.
%
%     Input:
%
%        n : order of the matrix;
%
%        d : diagonal elements of the input matrix;
%
%        e : subdiagonal elements of the input matrix
%            in its first n-1 positions.  e(n) is arbitrary.
%
%      Output:
%
%        d : eigenvalues in ascending order. if an error exit
%            is made, the eigenvalues are correct but unordered
%            for indices 1, 2, ..., ierr-1;
%
%        e : has been destroyed;
%
%        z : first components of the orthonormal eigenvectors
%            of the symmetric tridiagonal matrix.  if an error exit is
%            made, z contains the eigenvectors associated with the stored
%            eigenvalues;
%
%        ierr : error flag set to
%             0,  for normal return,
%             j,  if the j-th eigenvalue has not been
%                 determined after 30 iterations.
%
%     ------------------------------------------------------------------
%
%                                                       Restyle 2025-09-29.
%
z    = zeros(1,n); 
z(1) = 1;
ierr = 0;
if n ~= 1 
%
   e(n) = 0;
   for l = 1:n
       j = 0;
%     :::::::::: look for small sub-diagonal element ::::::::::
       ivar = 1;
       while ivar == 1
%
            for  m = l:n
                 if m == n
                    break
                 end
                 if abs(e(m)) <= eps*(abs(d(m)) + abs(d(m+1)))
                    break
                 end
            end 
%
            p = d(l);
            if  m ~= l
                if j == 30
                   ierr = l;
                   return
                end
                j = j + 1;
%     :::::::::: form shift ::::::::::
                g   = (d(l+1) - p) / (2*e(l));
                r   = sqrt(g*g+1);
                if g < 0
                   ggg = -abs(r);
                else
                   ggg = abs(r);
                end
                g   = d(m) - p + e(l) / (g + ggg);
                s   = 1;
                c   = 1;
                p   = 0;
                mml = m - l;
%     :::::::::: for i=m-1 step -1 until l do -- ::::::::::
               for ii = 1:mml  
                   i      = m - ii;
                   f      = s*e(i);
                   b      = c*e(i);
                   if abs(f) >= abs(g)
                      c      = g / f;
                      r      = sqrt(c*c+1);
                      e(i+1) = f*r;
                      s      = 1 / r;
                      c      = c*s;
                   else
                      s      = f / g;
                      r      = sqrt(s*s+1);
                      e(i+1) = g*r;
                      c      = 1 / r;
                      s      = s*c;
                   end 
                   g      = d(i+1) - p;
                   r      = (d(i) - g)*s + 2*c*b;
                   p      = s*r;
                   d(i+1) = g + p;
                   g      = c*r - b;
%     :::::::::: form first component of vector ::::::::::
                   f      = z(i+1);
                   z(i+1) = s*z(i) + c*f;
                   z(i) = c*z(i) - s*f;
               end 
               d(l) = d(l) - p;
               e(l) = g;
               e(m) = 0;
            else
               ivar = 2;
            end  
       end % while
   end  
%  :::::::::: order eigenvalues and eigenvectors ::::::::::
   for  ii = 2:n
        i = ii - 1;
        k = i;
        p = d(i);
        for  j = ii:n
             if d(j) <= p
                k = j;
                p = d(j);
             end 
        end  
        if k ~= i
           d(k) = d(i);
           d(i) = p;
           p    = z(i);
           z(i) = z(k);
           z(k) = p;
         end  
   end 
%     :::::::::: set error -- no convergence to an
%                eigenvalue after 30 iterations ::::::::::
end  
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = givens_GV(a,b)
%
%  G = givens_GV(a,b)
%
%  computation of the Givens rotation
%
%  G=[c -s; s c];
%
%  such that
%
%  G*[a;b] = [sqrt(a^2+b^2);0].
%
%  Input:
%
%    a : scalar;
%    b : scalar.
%
%  Output:
%    G : the Givens matrix.
%
%                                                       Restyle 2025-09-29.
%
if b==0
   c = 1; s = 0;
else
   if abs(b)>abs(a)
      tau = -a/b;
      s   = 1/sqrt(1+tau^2);
      c   = s*tau;
   else
      tau = -b/a;
      c   = 1/sqrt(1+tau^2);
      s   = c*tau;
   end
end
G = [c -s; s c];
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EOF %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%