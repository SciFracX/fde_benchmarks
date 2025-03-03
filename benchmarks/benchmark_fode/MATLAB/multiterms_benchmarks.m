clc

H=[2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7)];% 2^(-8)];% setting step-size

y0 = [1, 1, -1];
alpha = [3 2.5 2 1 0.5 0];
lambda = [1, 1, 1, 4, 1, 4];
f_fun = @(t, y) 6*cos(t);
J_fun = @(t, y) 0;
t0 = 0;
T = 100;
b=1;nb=0;

exa = @(t) sqrt(2).*sin(t + pi/4);

%benchmarking
Bench1=zeros(length(H),2);
Bench2=zeros(length(H),2);
Bench3=zeros(length(H),2);
Bench4=zeros(length(H),2);
Bench5=zeros(length(H),2);

for i=1:length(H)
    h=H(i);
%computing the time
t_series = t0:h:T;
u = f_fun(t_series, 0);
Bench1(i,1) = timeit(@() mt_fde_pi1_ex(alpha,lambda,f_fun,t0,T,y0,h));
Bench2(i,1) = timeit(@() mt_fde_pi1_im(alpha,lambda,f_fun,J_fun,t0,T,y0,h));
Bench3(i,1) = timeit(@() mt_fde_pi2_im(alpha,lambda,f_fun,J_fun,t0,T,y0,h));
Bench4(i,1) = timeit(@() mt_fde_pi12_pc(alpha,lambda,f_fun,t0,T,y0,h));
Bench5(i,1) = timeit(@() fode_caputo9(lambda,alpha,b,nb,y0,u,t_series));

%computing the error
[t1,y1]=mt_fde_pi1_ex(alpha,lambda,f_fun,t0,T,y0,h);
[t2,y2]=mt_fde_pi1_im(alpha,lambda,f_fun,J_fun,t0,T,y0,h);
[t3,y3]=mt_fde_pi2_im(alpha,lambda,f_fun,J_fun,t0,T,y0,h);
[t4,y4]=mt_fde_pi12_pc(alpha,lambda,f_fun,t0,T,y0,h);
y5 = fode_caputo9(lambda,alpha,b,nb,y0,u,t_series,2);

exact = exa(t1);
Bench1(i,2)=norm((y1-exact));
Bench2(i,2)=norm((y2-exact));
Bench3(i,2)=norm((y3-exact));
Bench4(i,2)=norm((y4-exact));
Bench5(i,2)=norm((y5'-exact));
end
%%
writematrix(Bench1,'MATLAB_MTPIEX.csv')
writematrix(Bench2,'MATLAB_MTPITrap.csv')
writematrix(Bench3,'MATLAB_MTPIRect.csv')
writematrix(Bench4,'MATLAB_MTPECE.csv')
writematrix(Bench5,'MATLAB_MTCAPUTO9.csv')