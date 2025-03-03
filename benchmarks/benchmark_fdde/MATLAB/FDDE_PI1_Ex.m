function [t, y] = FDDE_PI1_Ex(alpha,g,tau,t0,T,phi,h)

% This code solves a nonlinear fractional delay differential equation (FDDE) 
% with one constant delay tau > 0 in the form
%
% D^alpha y(t) = g(t, y(t), y(t-tau))
% y(t) = phi(t)  t0-tau <= t <= t0 
%
% where D^alpha is the fractional Caputo derivative of order 0 < alpha < 1.
%
% As initial data it must be provided not just a single value but a
% whole function phi(t) for t in the interval [t0-tau, t0]. 
%
% The FDDE is solved by an explicit rectangular Product-Integration (PI)
% scheme suitably modified to solve FDEs with one constant delay.
%
% Usage  [t, y] = FDDE_PI1_Ex(alpha,g,tau,t0,T,phi,h)
%
% Further information about this code are available in the Section 6 of the
% paper [1]; please, cite this code as [1]. 
%
% -------------------------------------------------------------------------
% Description of input parameters
% -------------------------------------------------------------------------
% alpha    : fractional order of the delay differential equation; it must
%            be 0 < alpha < 1  
% g        : function handle evaluating the vector field g(t,y(t),y(t-tau)) 
% tau      : constant delay; it must be tau > 0
% t0, T    : starting and final time of the interval of integration
% phi      : function handle for the initial data phi(t) 
% h        : integration step-size; it must be selected such that h <= tau
%
% -------------------------------------------------------------------------
% References and other information
% -------------------------------------------------------------------------
%
%  [1] Garrappa R., Kaslik E.: On initial conditions for fractional delay
%  differential equations, Communications in Nonlinear Science and
%  Numerical Simulation, 2020, 90, 105359, doi: 10.1016/j.cnsns.2020.105359  
%  Version of August, 10 2020
%
%  Please, report any problem or comment to :
%          roberto dot garrappa at uniba dot it
%
%  Copyright (c) 2020
%
%  Author:
%   Roberto Garrappa (University of Bari, Italy)
%   roberto dot garrappa at uniba dot it
%   Homepage: https://www.dm.uniba.it/members/garrappa
%
% -------------------------------------------------------------------------


if alpha <= 0 || alpha >= 1 || abs(imag(alpha)) > 0
   error('MATLAB:OrderOutsideRange',...
        ['The order ALPHA of the FDDE must 0 < ALPHA < 1. The value ' ...
        'ALPHA = %f is not supported.'], alpha);
end

if  tau < h
   error('MATLAB:StepSizeDelay',...
        ['The step-size H must be not greater than the constant delay ' ...
        'TAU. Reduce the step-size H such that H < %e.'], tau);
end   

N = ceil((T-t0)/h) ; 

% Evaluation of the grid points
t = t0 + h*(0:N) ;

% Evaluation of the PI1 coefficients 
nn_al = (0 : N).^alpha ;
b = [0, nn_al(2:end) - nn_al(1:end-1)]/gamma(alpha+1) ;
h_al = h^alpha ; 

% Evaluation of the solution by the explicit PI rule of order 1
y0 = phi(t0) ; 
y = zeros(length(y0),N+1) ; f = zeros(length(y0),N+1) ;
y(:,1) = y0 ;

for n = 1 : N
     
    % Evaluation of the vector field at the previous time tn-h
    tnm1 = t(n) ; 
    if tnm1 <= tau
        % Evaluation by means of the initial function
        y_nm1_tau = phi(tnm1-tau) ;
    else
        % Evaluation of the solution at the delyed time
        nm1_tau1 = floor((tnm1-tau)/h) ; nm1_tau2 = ceil((tnm1-tau)/h) ;
        if nm1_tau1 == nm1_tau2
            % The delayed solution belongs to the mesh: its value is used
            y_nm1_tau = y(nm1_tau1+1) ;
        else
            % The delayed solution does not belong to the mesh: its
            % interpolation in the interval between the closest mesh points
            % is used
            tt0 = t(nm1_tau1+1) ; tt1 = t(nm1_tau1+2) ;
            yy0 = y(nm1_tau1+1) ; yy1 = y(nm1_tau1+2) ; 
            y_nm1_tau = ((tnm1-tau)-tt0)/(tt1-tt0)*yy1 + ((tnm1-tau)-tt1)/(tt0-tt1)*yy0 ;
        end
    end
    f(:,n) = g(tnm1,y(:,n),y_nm1_tau) ;
        
    f_mem = 0 ;
    for j = 0 : n-1
        f_mem = f_mem + b(n-j+1)*f(:,j+1) ;
    end
    y(:,n+1) = y0 + h_al*f_mem ;
end

end





