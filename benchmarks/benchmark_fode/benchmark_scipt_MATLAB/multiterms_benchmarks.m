clc

% Step-size grid used for runtime and error benchmarking.
H=[2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7)];% 2^(-8)];% setting step-size

% Initial data for the multi-term equation.
% For this problem, y0 stores initial conditions needed by different orders.
y0 = [1, 1, -1];

% Orders of derivatives in the multi-term Caputo model.
alpha = [3 2.5 2 1 0.5 0];

% Coefficients multiplying each derivative term.
lambda = [1, 1, 1, 4, 1, 4];

% External forcing term and Jacobian (Jacobian is zero because f does not depend on y).
f_fun = @(t, y) 6*cos(t);
J_fun = @(t, y) 0;

% Time interval.
t0 = 0;
T = 100;

% Parameters for fode_caputo9 interface.
b=1;nb=0;

% Closed-form solution used as reference for error norms.
exa = @(t) sqrt(2).*sin(t + pi/4);

% Benchmark result matrices:
% column 1 -> runtime in seconds, column 2 -> error norm.
Bench1=zeros(length(H),2);
Bench2=zeros(length(H),2);
Bench3=zeros(length(H),2);
Bench4=zeros(length(H),2);
Bench5=zeros(length(H),2);

% Sweep over all step sizes.
for i=1:length(H)
    h=H(i);

% Build uniform grid and forcing samples for fode_caputo9.
t_series = t0:h:T;
u = f_fun(t_series, 0);

% Runtime benchmarking for each solver.
Bench1(i,1) = timeit(@() mt_fde_pi1_ex(alpha,lambda,f_fun,t0,T,y0,h));
Bench2(i,1) = timeit(@() mt_fde_pi1_im(alpha,lambda,f_fun,J_fun,t0,T,y0,h));
Bench3(i,1) = timeit(@() mt_fde_pi2_im(alpha,lambda,f_fun,J_fun,t0,T,y0,h));
Bench4(i,1) = timeit(@() mt_fde_pi12_pc(alpha,lambda,f_fun,t0,T,y0,h));
Bench5(i,1) = timeit(@() fode_caputo9(lambda,alpha,b,nb,y0,u,t_series));

% Compute numerical solutions for error evaluation at current h.
[t1,y1]=mt_fde_pi1_ex(alpha,lambda,f_fun,t0,T,y0,h);
[t2,y2]=mt_fde_pi1_im(alpha,lambda,f_fun,J_fun,t0,T,y0,h);
[t3,y3]=mt_fde_pi2_im(alpha,lambda,f_fun,J_fun,t0,T,y0,h);
[t4,y4]=mt_fde_pi12_pc(alpha,lambda,f_fun,t0,T,y0,h);
y5 = fode_caputo9(lambda,alpha,b,nb,y0,u,t_series,2);

% Reference trajectory on the time grid returned by method 1.
exact = exa(t1);

% Error norms against the analytic solution.
Bench1(i,2)=norm((y1-exact));
Bench2(i,2)=norm((y2-exact));
Bench3(i,2)=norm((y3-exact));
Bench4(i,2)=norm((y4-exact));
Bench5(i,2)=norm((y5'-exact));
end
%%
% Export benchmark tables to CSV.
% Each CSV stores [runtime, error] for all step sizes in H.
writematrix(Bench1,'MATLAB_MTPIEX.csv')
writematrix(Bench2,'MATLAB_MTPITrap.csv')
writematrix(Bench3,'MATLAB_MTPIRect.csv')
writematrix(Bench4,'MATLAB_MTPECE.csv')
writematrix(Bench5,'MATLAB_MTCAPUTO9.csv')