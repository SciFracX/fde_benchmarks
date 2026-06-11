clc

% Step-size grid for runtime/accuracy benchmarking.
H=[2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7)];% 2^(-8)];% setting step-size

% Fractional order for this linear single-term test problem.
alpha = 0.8;

% Linear right-hand side D^alpha y = -10*y.
f_fun = @(t,y) -10*y;

% Jacobian with respect to y (constant for this linear model).
J_fun = @(t,y) -10 ;

% Time interval and initial condition.
t0 = 0.0;
T = 5.0;
y0 = 1.0;

% Exact solution using Mittag-Leffler function.
exa = @(t) ml(-10*t.^alpha, alpha);

% Scalar wrapper required by nlfode_vec interface.
% k is unused because the model has only one state.
function y=fun(t,x,k)
        y = -10*x(1) ;
end

% Benchmark result matrices:
% column 1 -> runtime in seconds, column 2 -> solution error norm.
Bench1=zeros(length(H),2);
Bench2=zeros(length(H),2);
Bench3=zeros(length(H),2);
Bench4=zeros(length(H),2);
Bench5=zeros(length(H),2);
Bench6=zeros(length(H),2);
Bench7=zeros(length(H),2);
Bench8=zeros(length(H),2);

% Sweep all step sizes and collect timing + error statistics.
for i=1:length(H)
    h=H(i);

% Runtime measurements for each solver.
Bench1(i,1) = timeit(@() fde_pi1_ex(alpha,f_fun,t0,T,y0,h));
Bench2(i,1) = timeit(@() fde_pi12_pc(alpha,f_fun,t0,T,y0,h));
Bench3(i,1) = timeit(@() fde_pi1_im(alpha,f_fun,J_fun,t0,T,y0,h));
Bench4(i,1) = timeit(@() fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h));
Bench5(i,1) = timeit(@() nlfode_vec(@fun,alpha,y0,h,T));
Bench6(i,1) = timeit(@() flmm2(alpha,f_fun,J_fun,t0,T,y0,h,[],1));
Bench7(i,1) = timeit(@() flmm2(alpha,f_fun,J_fun,t0,T,y0,h,[],2));
Bench8(i,1) = timeit(@() flmm2(alpha,f_fun,J_fun,t0,T,y0,h,[],3));

% Numerical solves at current h for error computation.
[t1,y1]=fde_pi1_ex(alpha,f_fun,t0,T,y0,h);
[t2,y2]=fde_pi12_pc(alpha,f_fun,t0,T,y0,h);
[t3,y3]=fde_pi1_im(alpha,f_fun,J_fun,t0,T,y0,h);
[t4,y4]=fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h);
[y5,t5]=nlfode_vec(@fun,alpha,y0,h,T);
[t6,y6]=flmm2(alpha,f_fun,J_fun,t0,T,y0,h,[],1);
[t7,y7]=flmm2(alpha,f_fun,J_fun,t0,T,y0,h,[],2);
[t8,y8]=flmm2(alpha,f_fun,J_fun,t0,T,y0,h,[],3);

% Reference solution sampled on the time grid returned by method 1.
exact = exa(t1);

% Error norms against the exact solution trajectory.
Bench1(i,2)=norm((y1-exact));
Bench2(i,2)=norm((y2-exact));
Bench3(i,2)=norm((y3-exact));
Bench4(i,2)=norm((y4-exact));
Bench5(i,2)=norm((y5'-exact));
Bench6(i,2)=norm((y6-exact));
Bench7(i,2)=norm((y7-exact));
Bench8(i,2)=norm((y8-exact));
end
%%
% Export benchmark tables to CSV.
% Each file stores two columns: [runtime, error].
writematrix(Bench1,'Linear_Singleterm_MATLAB_PIEX.csv')
writematrix(Bench2,'Linear_Singleterm_MATLAB_PECE.csv')
writematrix(Bench3,'Linear_Singleterm_MATLAB_PIRect.csv')
writematrix(Bench4,'Linear_Singleterm_MATLAB_PITrap.csv')
writematrix(Bench5,'Linear_Singleterm_MATLAB_NLFODE_VEC.csv')
writematrix(Bench6,'Linear_Singleterm_MATLAB_Trapzoid.csv')
writematrix(Bench7,'Linear_Singleterm_MATLAB_NewtonGregory.csv')
writematrix(Bench8,'Linear_Singleterm_MATLAB_BDF.csv')