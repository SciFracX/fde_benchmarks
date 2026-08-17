clc

H=[2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7) 2^(-8)];% setting step-size

alpha = [0.5, 0.2, 0.6];
f_fun = @(t,y) [(((y(2)-0.5).*(y(3)-0.3)).^(1/6) + sqrt(t))/sqrt(pi);
        gamma(2.2)*(y(1)-1);
        gamma(2.8)/gamma(2.2)*(y(2)-0.5)];
J_fun = @(t,y) [0, (y(2)-0.5).^(-5/6).*(y(3)-0.3).^(1/6)/6/sqrt(pi), (y(2)-0.5).^(1/6).*(y(3)-0.3).^(-5/6)/6/sqrt(pi);
        gamma(2.2), 0 , 0;
        0 , gamma(2.8)/gamma(2.2) , 0];
t0 = 0;
T = 5;
y0 = [ 1 ; 0.5 ; 0.3 ] ;
exa = @(t) [t + 1; t.^1.2 + 0.5; t.^1.8 + 0.3];

%benchmarking
Bench1=zeros(length(H),2);
Bench2=zeros(length(H),2);
Bench3=zeros(length(H),2);
Bench4=zeros(length(H),2);
Bench5=zeros(length(H),2);

for i=1:length(H)
    h=H(i);
    N=round((T-t0)/h);
%computing the time
Bench1(i,1) = timeit(@() fde_pi1_ex(alpha,f_fun,t0,T,y0,h));
Bench2(i,1) = timeit(@() fde_pi12_pc(alpha,f_fun,t0,T,y0,h));
Bench3(i,1) = timeit(@() fde_pi1_im(alpha,f_fun,J_fun,t0,T,y0,h));
Bench4(i,1) = timeit(@() fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h));
Bench5(i,1) = timeit(@() nlfode_vec(@fun,alpha,y0,h,T));
%computing the error
[t1,y1]=fde_pi1_ex(alpha,f_fun,t0,T,y0,h);
[t2,y2]=fde_pi12_pc(alpha,f_fun,t0,T,y0,h);
[t3,y3]=fde_pi1_im(alpha,f_fun,J_fun,t0,T,y0,h);
[t4,y4]=fde_pi2_im(alpha,f_fun,J_fun,t0,T,y0,h);
[y5,t5]=nlfode_vec(@fun,alpha,y0,h,T);

exact = exa(t3);
Bench1(i,2)=norm((y1-exact), Inf);
Bench2(i,2)=norm((y2-exact), Inf);
Bench3(i,2)=norm((y3-exact), Inf);
Bench4(i,2)=norm((y4-exact), Inf);
Bench5(i,2)=norm((y5'-exact), Inf);
end
%%
data_dir = '/Users/quqingyu/SciFracX/paper/benchmarks/data';
if ~isfolder(data_dir)
    mkdir(data_dir)
end
% Convert benchmark matrices to tables with column names
T1 = array2table(Bench1, 'VariableNames', {'time', 'error'});
T2 = array2table(Bench2, 'VariableNames', {'time', 'error'});
T3 = array2table(Bench3, 'VariableNames', {'time', 'error'});
T4 = array2table(Bench4, 'VariableNames', {'time', 'error'});
T5 = array2table(Bench5, 'VariableNames', {'time', 'error'});

% Write CSV files with headers
writetable(T1, fullfile(data_dir, 'MATLAB_PIEX.csv'));
writetable(T2, fullfile(data_dir, 'MATLAB_PECE.csv'));
writetable(T3, fullfile(data_dir, 'MATLAB_PIRect.csv'));
writetable(T4, fullfile(data_dir, 'MATLAB_PITrap.csv'));
writetable(T5, fullfile(data_dir, 'MATLAB_NLFODE_VEC.csv'));

% Component-wise wrapper required by the nlfode_vec interface.
function value = fun(t,x,k)
    switch k
        case 1
            value = (((x(2)-0.5).*(x(3)-0.3)).^(1/6) ...
                + sqrt(t))/sqrt(pi);
        case 2
            value = gamma(2.2)*(x(1)-1);
        case 3
            value = gamma(2.8)/gamma(2.2)*(x(2)-0.5);
    end
end

% Unified callback for the 13-state commensurate formulation used by
% fhbvm and fhbvm2.  The common fractional order is 0.1.
function value = g_fun(t,y,~)
    if nargin == 0
        value = 0.1;
    elseif nargin == 2
        value = zeros(size(y));
        value(:,1:4) = y(:,2:5);
        product = (y(:,6)-0.5).*(y(:,8)-0.3);
        value(:,5) = (abs(product).^(1/6)+sqrt(t(:)))/sqrt(pi);
        value(:,6) = y(:,7);
        value(:,7) = gamma(2.2).*(y(:,1)-1);
        value(:,8:12) = y(:,9:13);
        value(:,13) = gamma(2.8)/gamma(2.2).*(y(:,6)-0.5);
    elseif nargin == 3
        value = zeros(13);
        value(1:4,2:5) = eye(4);
        value(6,7) = 1;
        value(7,1) = gamma(2.2);
        value(8:12,9:13) = eye(5);
        value(13,6) = gamma(2.8)/gamma(2.2);

        % The two omitted entries of the exact Jacobian are singular at
        % t=0.  Zero is used as a regularized simplified-Newton Jacobian.
        value(5,6) = 0;
        value(5,8) = 0;
    else
        error('g_fun:InvalidInputCount', ...
            'g_fun expects 0, 2, or 3 input arguments.');
    end
end
