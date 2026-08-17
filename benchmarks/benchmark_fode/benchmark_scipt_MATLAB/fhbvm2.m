function [t,y,etim] = fhbvm2( fun, y0, T, N, n, nu, txton )
%--------------------------------------------------------------------------
%
%   [t,y,etim] = fhbvm2( fun, y0, T, N, n, nu )
%
%   Solves the initial value problem for fractional differential equations 
%   of Caputo type:
%
%    y^(alpha) = f(t,y),   0<=t<=T,   
%
%    y^(i)(0)  = y_i \in R^m, i = 0,...,[alpha]-1,
%
%   with [alpha] = ceil(alpha).  (REMARK: alpha>0 and not an integer)
%
%    Input:
%    ======
%    - y0    : matrix  [alpha] x m  with the initial conditions. 
%              In more detail,
%
%                     y0(i+1,:) = (y_i).',  i = 0,...,[alpha]-1.
%
%    - T     : final integration time.
%
%    - fun   : a string containing a function identifier, 
%              or a function_handle, such that:
%
%              alpha   = fun()
%
%              f(t,y)  = fun( t, y ) (REMARK: MUST WORK  IN  VECTOR MODE,
%                                             i.e., if [ell,m] = size(y),
%                                             and  length(t) == ell, then
%                                             fun(t,y)   must  return  an  
%                                             ell x m  matrix)
%
%              f'(t,y) = fun( t, y, 1 ) (i.e., computes the Jacobian of f).
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
%    [5] Website: https://people.dimai.unifi.it/brugnano/fhbvm/
%
%
%                                                          Rel. 2025-07-28.
%
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% flag for the final output
data.txt  = nargin==7; 
txt       = data.txt;
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
data.fun  = fun;
data.k    = 22;   % fixed at this moment
data.s    = 22;   % fixed at this moment
k         = data.k;
s         = data.s;
data.alfa = feval(fun);
alfa      = data.alfa;
data.y0   = y0;
data.m    = size(y0,2);
%m        = data.m;  %%% not needed here
data.T    = T;       %%% width of the integration interval
data.n    = n;
data.n1   = n1;
if n1==1, r = 2; else, r = n1/(n1-1); end
data.r    = r;    
data.h    = T/n;
h         = data.h;  %%% timestep in the uniform mesh
numin     = ceil(log(11)/log(r));  %%%%%%% guarantees hnu <= 1.1*h %%%%%%%%
nup       = numin>nu && nu>n1;
if nup
    nu    = numin; 
end    
data.nu   = nu;
data.h0   = n1*h*(r-1)/(r^nu-1);   %%% initial timestep in the graded mesh
h0        = data.h0;
hnu       = r^(nu-1)*h0;           %%% last timestep in the graded mesh
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
data.Jjc  = [];
data.NaN  = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initialize execution time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% compute the Jacobi abscissae, weights, etc.
[x,beta,PO,Is,~,Xs,abd] = jacobidat( k, s, alfa );
x         = x(:);
x1        = [x; 1];          % the abscissa 1 is needed to advance the step
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
kGauss    = 30;   % order of the Gauss-Legendre quadrature used == 2*kGauss
[Ps,glc,glb] = makePs( abd, kGauss ); % Gauss-Legendreabscissae and weights
data.kGauss  = kGauss;                %     and corresponding Jacobi matrix
data.Ps      = Ps;     
data.glc     = glc;    
data.glb     = glb;  
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
if nu>1
   Jj = zeros(s,(k+1)*(nu-1));  %%% graded mesh
else
   Jj = [];
end   
ind   = 0;
for j = 1:nu-1  %% graded mesh
    for i = 1:k+1
        ind       = ind+1;
        M         = (r^j-1)/(r-1) + x1(i)*r^j;
        Jj(:,ind) = Jjalfa1( M, abd, x, beta, alfa, Ps, glc, glb );       
    end    
end
data.Jj  = Jj;
if n-n1>1
   Jjc = zeros(s,(k+1)*(n-n1)); %%% uniform mesh mesh
else
   Jjc = [];
end   
ind   = 0;
for j = 1:n-n1    %% (r^nu-1)/(r-1) h1 = n1*h
    for i = 1:k+1
        ind        = ind+1;
        M          = j + x1(i);
        Jjc(:,ind) = Jjalfa1( M, abd, x, beta, alfa, Ps, glc, glb );     
    end    
end
data.Jjc = Jjc;
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
   fprintf('Solver: FHBVM(%i,%i) \n',k,s)
   disp(linea)
   if nup
      fprintf('nu updated to %i \n',nu)
   end    
   if nu>n1
      fprintf('initial graded mesh with %i points \n',nu)
      fprintf('     initial timestep = %10.2e \n',h0)
      fprintf('     final   timestep = %10.2e \n',hnu)
      fprintf('     and parameter  r = %g \n',r)
      fprintf('\n')
   end
   if n>n1
      fprintf('%i uniform points with timestep h = %10.2e \n',n-n1+uno,h)
   end   
   disp(linea)
   fprintf( ' %i blended steps and %i fixed-point steps\n', ...
                                              blend, length(t)-blend-1 ); 
   fprintf( '    total iterations solver: \n' );
   fprintf( '                      blended = %i \n', it(1) );
   fprintf( '                  fixed-point = %i \n', it(2) );
   if flag_NaN
      fprintf(' NaN occurred at t = %g \n', t(end) );
   else
      fprintf(' execution regularly ended at t = %g \n',t(end)) 
   end   
   disp(linea)
   fprintf( ' execution time = %g sec \n', etim );
   disp(linea)
elseif flag_NaN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% an error message is issued
      fprintf(' NaN occurred at t = %g \n', t(end) );  
end
%-------------------------------------------------------------------------%
return

function [t,y,it,blend,flag_NaN] = fhbvmdata( data )  
%
%  [t,y,it,blend,flag_NaN] = fhbvmdata( data );
%
%  Solves the system of fractional ODEs, for alfa>0, 
%
%         D^alfa y(t) = f(t,y), t \in [0,T], 
%            y^(i)(0) = y_0^i \in R^m, i = 0,...,[alfa]-1,         
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
%                        alfa    = fun ; 
%
%                        f(t,y)  = fun(t,y)   (has to work in vector mode);
%
%                        f'(t,y) = fun(t,y,1)  (vector mode not needed);
%
%            data.k    - parameter k of FHBVM(k,s);
%            data.s    - parameter s of FHBVM(k,s);
%            data.alfa - order alfa of the derivative; 
%            data.y0   - initial condition(s);
%            data.m    - dimension of the problem;
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
%            data.Jj   - auxiliary integrals for the memory terms on the 
%                        graded mesh; 
%            data.Jjc  - auxiliary integrals for the memory terms on the
%                        uniform mesh;
%            data.Ps   - Jacobi polynomials evaluated at the GL abscissae;
%            data.glc  - GL abscissae of enough high order;
%            data.glb  - GL weights; 
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
%   flag_NaN - flag set to 1 in case NaN occur.
%
%--------------------------------------------------------------------------
%                                                          Rel. 2025-04-07.
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fun      = data.fun;                 %
y0       = data.y0;                  %
k        = data.k;                   %
s        = data.s;                   %
alfa     = data.alfa;                %
T        = data.T;                   %
n        = data.n;                   %
n1       = data.n1;                  %
nu       = data.nu;                  %
h0       = data.h0;                  %
r        = data.r;                   %
h        = data.h;                   %
x        = data.x;                   %
x1       = data.x1;                  %
PO       = data.PO;                  %
Is       = data.Is;                  %
fattb    = data.fatt;                %
gam      = data.gam;                 %
CC       = data.CC;                  %
Jj       = data.Jj;                  %
Jjc      = data.Jjc;                 %
m        = data.m; % == size(y0,2);  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hh       = [h0*cumprod( [1; r*ones(nu-1,1)] ); h*ones(n-n1,1)]; % timesteps
t        = [0; cumsum(hh)];                                     % and  mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if abs(T-t(end))>max(100,2*n)*(1+T)*eps    %
   disp('-------------------------------') %
   disp('fhbvmdata: mesh not compatible.') %  
   disp('-------------------------------') %  
   error('fhbvmdata: exit')                %    mesh consistency check
else                                       %    and correction  
   cT = T/t(end);                          %
   h  = h*cT;                              %
   t  = t*cT;                              %
end                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y        = zeros(nu+n-n1+1,m);   % y(1,:)        initial condition 
y(1,:)   = y0(1,:);              % y(2:nu+1,:)   graded mesh 
fatt     = 1/gamma(alfa);        % y(nu+2:end,:) uniform mesh
fatt1    = fatt/alfa;
halfa    = hh.^alfa;   % h_i^alfa on the graded mesh
half     = h^alfa;     % h_i^alfa on the uniform mesh
it       = zeros(1,2);
%
%  gamJ(:,(i-1)*s:i*s) = Jacobi-Fourier coefficients at step i, graded mesh
%  gamJc = Jacobi-Fourier coefficients at step i, uniform mesh
%                          
gamJ     = zeros(m,nu*s);         % on the graded mesh  
gamJc    = zeros(m,(n-n1)*s);     % on the uniform mesh
%-------------------------------------------------------------------------%
Im       = eye(m);
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
    blendflag       = L*fattb*halfa(i) > small;
    if blendflag
       teta         = inv( Im -halfa(i)*gam*J0.' );    % transpose Jacobian
       blend        = blend + 1;
    else
       teta         = [];
    end   
    tcor            = t(i) + hh(i)*x;
    tcor1           = [tcor; t(i+1)];
    fi0             = evalfi( y0, tcor1, x1, Jj(:,1:indk), gamJ(:,1:inds));   
    % fi0 = evalfi( y0, tcor1, x1, Jj(:,1:indk), gamJ); % do not improve %%   
    [yi,gami,iti,flag_NaN] = step( fun, tcor, fi0, halfa(i)*Is, PO, ...
                                             fatt1*halfa(i), CC, teta, L );
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
if flag_NaN, return, end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NaN: exit
inds  = 0;
indk  = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uniform mesh
for j = n1+uno:n  
    i = i+1;       
    J0              = feval( fun, t(i), y(i,:), 1 ); 
    L               = min( norm(J0,1), norm(J0,inf) );
    blendflag       = L*fattb*half > small;
    if blendflag
       teta         = inv( Im -half*gam*J0.' );    % transpose Jacobian
       blend        = blend + 1;
    else
       teta         = [];
    end   
    tcor            = t(i) + h*x;
    tcor1           = [tcor; t(i+1)];
    fi0             = evalfi(y0,tcor1, x1, Jjc(:,1:indk), gamJc(:,1:inds));
    % fi0 = evalfi(y0,tcor1, x1, Jjc(:,1:indk), gamJc); % do not improve %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         sum of the contribution of the graded mesh to the memory term
    fi0             = fi0 + evalJjinu( j, gamJ, data );
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [yi,gami,iti,flag_NaN] = step( fun, tcor, fi0, half*Is, PO, ...
                                                 fatt1*half, CC, teta, L );
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
    y(i+1,:)         = yi;
    gamJc(:,inds+s1) = ( fatt*half )*gami.';
    inds             = inds + s;
    indk             = indk + k+1;
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
%                                                          Rel. 2025-06-07.
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
kk1   = 1:k1;
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

function [y,gam,it,fNaN] = step( fun, t, fi0, Is, PO, ghalfa, CC, teta, L )
%
%   Integration step for the local problem.
%                                                          Rel. 2024-03-09.
%
[k1,m] = size(fi0);
k      = k1-1;
Y0     = fi0(1:k,:);
tol0   = 2*eps*(1+L);
tol    = 100*tol0;
itmax  = 1000*(m+1)*k1;
y      = Y0;
scal   = 1./( 1 + abs(Y0) );
err    = Inf;
gam    = 0;   % Fourier coefficients
noblen = isempty(teta); 
for it = 1:itmax
    yold = y;
    f    = feval(fun,t,y);
    if noblen % fixed-point iteration
       gam  = PO*f;
    else      % blended iteration
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
   y = fi0(k1,:) + ghalfa*gam(1,:); % new approximation
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
%       - k,s  : parameters of the FHBVM(k,s) method, 1<=s<=k;
%       - alfa : index of the fractional derivative (alfa>0 and not in IN);
%  OUTPUT:
%       - x    : k Gauss-Jacobi abscissae in [0,1], relative to the 
%                weight function
%                               w(x) = alfa*(1-x)^alfa;
%       - beta : corresponding quadrature weights; % sum(beta)=1 %
%       - PO   : PO = P_s^T * Omega, with P_s(i,j) = P_{j-1}(x_i)
%                                    and Omega = diag(beta);
%       - Is   : matrix I_s(i,j) = I^alfa P_{j-1}(x_i); % i = 1,...,k,
%       - Ps   : matrix P_s(i,j) = P_{j-1}(x_i);        % j = 1,...,s.
%       - Xs   : = PO*Is;
%       - abd  : coefficients [a b d] for the recurrence:
%               P_{j+1}(c)=(a(j)*c-b(j))*P_j(c)-d(j)*P_{j-1}(c), j=0,..s-2,
%               P_0(c)==1, P_{-1}(c)==0.
%
%                                                          Rel. 2025-03-09.
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
    warning('jacobidat: are abscissae accurate? please check') % keyboard 
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

function [Ps,x,beta] = makePs( abd, k )
%
%   Compute the k Gauss-Legendre nodes (x) and weights (beta), along with
%   the Jacobi matrix evaluated at the Gauss-abscissae (Ps).
%   abd contains the coefficients of the Jacobi recursion.
%
%                                                          Rel. 2025-02-07.
%   REMARK: order of the Gauss-Legendre quadrature = 2*k
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
%---------------------------------------------------------------------- EOF