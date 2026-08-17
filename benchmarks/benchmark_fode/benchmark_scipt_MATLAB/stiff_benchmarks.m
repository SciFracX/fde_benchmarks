clc

% Step-size grid used in the convergence and performance sweep.
H=[2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7) 2^(-8)];% setting step-size

% Fractional orders for the 3-component system.
alpha = [0.5, 0.5, 0.5];

% Scalar order used in some single-order formulations.
beta = 0.5;

% Coefficients of the manufactured polynomial-in-time exact solution.
a1 = 0.5; a2 = 0.8; a3 = 1; a4 = 1; a5 = 1; a6 = 1;

% Exponents used by the manufactured exact solution.
sigma = [beta, 2*beta, 1+beta, 5*beta, 2, 2+beta];

% Linear stiff dynamics matrix blocks.
A = [-10000 0 1;
     -0.05 -0.08 -0.2;
     1 0 -1];
B = [-0.6 0 0.2;
     -0.1 -0.2 0;
     0 -0.5 -0.8];

% Initial condition at t = 0.
u0 = [1; 1; 1];

% Gamma-ratio helper appearing in Caputo derivative of power functions.
big_gamma = @(k) gamma(sigma(k)+1)/gamma(sigma(k)-beta+1);

% Manufactured forcing term g(t) so the selected exact solution is known.
g = @(t) [a1*big_gamma(1)*t^(sigma(1)-beta) + a2*big_gamma(2)*t^(sigma(2)-beta); a3*big_gamma(3)*t^(sigma(3)-beta) + a4*big_gamma(4)*t^(sigma(4)-beta); a5*big_gamma(5)*t^(sigma(5)-beta) + a6*big_gamma(6)*t^(sigma(6)-beta)] - (A + B) * [a1*t^sigma(1) + a2*t^sigma(2) + u0(1); a3*t^sigma(3) + a4*t^sigma(4) + u0(2); a5*t^sigma(5) + a6*t^sigma(6) + u0(3)];

% System right-hand side used by most benchmarked solvers.
f_fun = @(t, y) A*y + B*y + g(t);

% Constant Jacobian for implicit methods that support user-supplied J.
J_fun = @(t, y) A + B;

% Time interval.
t0 = 0;
T = 1;

% Analytic solution used to compute terminal errors.
exa = @(t) [a1*t.^sigma(1) + a2*t.^sigma(2) + u0(1); a3*t.^sigma(3) + a4*t.^sigma(4) + u0(2); a5*t.^sigma(5) + a6*t.^sigma(6) + u0(3)];

% fhbvm2 uses ten graded points over the first interval of width h and a
% uniform mesh with stepsize h over the remainder of the integration span.
fhbvm2_graded_span = 1;
fhbvm2_graded_points = 10;

% Benchmark result arrays: column 1 = runtime (seconds), column 2 = error.
Bench1=zeros(length(H),2);
Bench2=zeros(length(H),2);
Bench3=zeros(length(H),2);
Bench4=zeros(length(H),2);
Bench5=zeros(length(H),2);
Bench6=zeros(length(H),2);
Bench7=zeros(length(H),2);
Bench8=zeros(length(H),2);
Bench9=zeros(length(H),2);

% Sweep all step sizes and benchmark each solver.
for i=1:length(H)
    h=H(i);
    N=round((T-t0)/h);

% Runtime measurement using timeit for fair timing.
Bench1(i,1) = timeit(@() fde_pi1_ex(alpha,f_fun,t0,T,u0,h));
Bench2(i,1) = timeit(@() fde_pi12_pc(alpha,f_fun,t0,T,u0,h));
Bench3(i,1) = timeit(@() fde_pi1_im(alpha,f_fun,J_fun,t0,T,u0,h));
Bench4(i,1) = timeit(@() fde_pi2_im(alpha,f_fun,J_fun,t0,T,u0,h));
Bench5(i,1) = timeit(@() flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],1));
Bench6(i,1) = timeit(@() flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],2));
Bench7(i,1) = timeit(@() flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],3));
Bench8(i,1) = timeit(@() nlfode_vec(@fun,alpha,u0,h,T));
Bench9(i,1) = timeit(@() pepc_nlfode(f_fun,beta,u0,h,T,1e-6,2));

% Solve once with each method at the current step size.
[t1,y1]=fde_pi1_ex(alpha,f_fun,t0,T,u0,h);
[t2,y2]=fde_pi12_pc(alpha,f_fun,t0,T,u0,h);
[t3,y3]=fde_pi1_im(alpha,f_fun,J_fun,t0,T,u0,h);
[t4,y4]=fde_pi2_im(alpha,f_fun,J_fun,t0,T,u0,h);
[t5,y5]=flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],1);
[t6,y6]=flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],2);
[t7,y7]=flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],3);
[y8,t8]=nlfode_vec(@fun,alpha,u0,h,T);
[y9,t9]=pepc_nlfode(f_fun,beta,u0,h,T,1e-6,2);

% Reference solution sampled at t-grid returned by method 1.
exact = exa(t1);

% Terminal/global error proxies (vector norms versus analytic trajectory).
Bench1(i,2)=norm((y1-exact), Inf);
Bench2(i,2)=norm((y2-exact), Inf);
Bench3(i,2)=norm((y3-exact), Inf);
Bench4(i,2)=norm((y4-exact), Inf);
Bench5(i,2)=norm((y5-exact), Inf);
Bench6(i,2)=norm((y6-exact), Inf);
Bench7(i,2)=norm((y7-exact), Inf);
Bench8(i,2)=norm((y8'-exact), Inf);
Bench9(i,2)=norm((y9'-exact), Inf);
end
%%
% Export per-solver benchmark tables.
% Each CSV stores [runtime, error] per step size in H.
data_dir = '/Users/quqingyu/SciFracX/paper/benchmarks/data';
if ~isfolder(data_dir)
    mkdir(data_dir)
end
% Export benchmark results with CSV headers: time,error

writetable(array2table(Bench1, 'VariableNames', {'time','error'}), fullfile(data_dir, 'Stiff_MATLAB_PIEX.csv'));
writetable(array2table(Bench2, 'VariableNames', {'time','error'}), fullfile(data_dir, 'Stiff_MATLAB_PECE.csv'));
writetable(array2table(Bench3, 'VariableNames', {'time','error'}), fullfile(data_dir, 'Stiff_MATLAB_PIRect.csv'));
writetable(array2table(Bench4, 'VariableNames', {'time','error'}), fullfile(data_dir, 'Stiff_MATLAB_PITrap.csv'));
writetable(array2table(Bench5, 'VariableNames', {'time','error'}), fullfile(data_dir, 'Stiff_MATLAB_Trapzoid.csv'));
writetable(array2table(Bench6, 'VariableNames', {'time','error'}), fullfile(data_dir, 'Stiff_MATLAB_NewtonGregory.csv'));
writetable(array2table(Bench7, 'VariableNames', {'time','error'}), fullfile(data_dir, 'Stiff_MATLAB_BDF.csv'));
writetable(array2table(Bench8, 'VariableNames', {'time','error'}), fullfile(data_dir, 'Stiff_MATLAB_NLFODE_VEC.csv'));
writetable(array2table(Bench9, 'VariableNames', {'time','error'}), fullfile(data_dir, 'Stiff_MATLAB_PEPC_NLFODE.csv'));

% Component-wise wrapper required by the nlfode_vec interface.
function value = fun(t,x,k)
    beta = 0.5;
    u0 = [1;1;1];
    sigma = [beta,2*beta,1+beta,5*beta,2,2+beta];
    coeff = [0.5,0.8,1,1,1,1];
    big_gamma = @(j) gamma(sigma(j)+1)/gamma(sigma(j)-beta+1);
    A = [-10000 0 1; -0.05 -0.08 -0.2; 1 0 -1];
    B = [-0.6 0 0.2; -0.1 -0.2 0; 0 -0.5 -0.8];
    exact = [coeff(1)*t^sigma(1)+coeff(2)*t^sigma(2)+u0(1); ...
        coeff(3)*t^sigma(3)+coeff(4)*t^sigma(4)+u0(2); ...
        coeff(5)*t^sigma(5)+coeff(6)*t^sigma(6)+u0(3)];
    derivative = [coeff(1)*big_gamma(1)*t^(sigma(1)-beta) ...
            + coeff(2)*big_gamma(2)*t^(sigma(2)-beta); ...
        coeff(3)*big_gamma(3)*t^(sigma(3)-beta) ...
            + coeff(4)*big_gamma(4)*t^(sigma(4)-beta); ...
        coeff(5)*big_gamma(5)*t^(sigma(5)-beta) ...
            + coeff(6)*big_gamma(6)*t^(sigma(6)-beta)];
    rhs = (A+B)*x+derivative-(A+B)*exact;
    value = rhs(k);
end

% Vectorized callback required by fhbvm and fhbvm2.
function value = g_fun(t,y,~)
    if nargin == 0
        value = 0.5;
    elseif nargin == 2
        A = [-10000 0 1; -0.05 -0.08 -0.2; 1 0 -1];
        B = [-0.6 0 0.2; -0.1 -0.2 0; 0 -0.5 -0.8];
        sigma = [0.5,1,1.5,2.5,2,2.5];
        coeff = [0.5,0.8,1,1,1,1];
        tt = t(:);
        exact = [coeff(1)*tt.^sigma(1)+coeff(2)*tt.^sigma(2)+1, ...
            coeff(3)*tt.^sigma(3)+coeff(4)*tt.^sigma(4)+1, ...
            coeff(5)*tt.^sigma(5)+coeff(6)*tt.^sigma(6)+1];
        derivative = [ ...
            coeff(1)*gamma(sigma(1)+1)/gamma(sigma(1)+0.5).*tt.^(sigma(1)-0.5) ...
                + coeff(2)*gamma(sigma(2)+1)/gamma(sigma(2)+0.5).*tt.^(sigma(2)-0.5), ...
            coeff(3)*gamma(sigma(3)+1)/gamma(sigma(3)+0.5).*tt.^(sigma(3)-0.5) ...
                + coeff(4)*gamma(sigma(4)+1)/gamma(sigma(4)+0.5).*tt.^(sigma(4)-0.5), ...
            coeff(5)*gamma(sigma(5)+1)/gamma(sigma(5)+0.5).*tt.^(sigma(5)-0.5) ...
                + coeff(6)*gamma(sigma(6)+1)/gamma(sigma(6)+0.5).*tt.^(sigma(6)-0.5)];
        forcing = derivative-exact*(A+B).';
        value = y*(A+B).'+forcing;
    elseif nargin == 3
        value = [-10000.6 0 1.2; -0.15 -0.28 -0.2; 1 -0.5 -1.8];
    else
        error('g_fun:InvalidInputCount', ...
            'g_fun expects 0, 2, or 3 input arguments.');
    end
end
