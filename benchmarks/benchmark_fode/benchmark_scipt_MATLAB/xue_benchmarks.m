clc

H=[2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7)];% 2^(-8)];% setting step-size

alpha = [0.5, 0.2, 0.6];
f_fun = @(t,y) [(((y(2)-0.5).*(y(3)-0.3)).^(1/6) + sqrt(t))/sqrt(pi);
        gamma(2.2)*(y(1)-1);
        gamma(2.8)/gamma(2.2)*(y(2)-0.5)];
J_fun = @(t,y) [0, (y(2)-0.5).^(-5/6).*(y(3)-0.3).^(1/6)/6/sqrt(pi), (y(2)-0.5).^(1/6).*(y(3)-0.3).^(-5/6)/6/sqrt(pi);
        gamma(2.2), 0 , 0;
        0 , gamma(2.8)/gamma(2.2) , 0];
t0 = 0;
T = 5;
y0 = [ 1 ; 0.500000001 ; 0.300000001 ] ;
exa = @(t) [t + 1; t.^1.2 + 0.5; t.^1.8 + 0.3];
function y=fun(t,x,k)
switch k
    case 1
        y=(((x(2)-0.5).*(x(3)-0.3)).^(1/6) + sqrt(t))/sqrt(pi);
    case 2
        y=gamma(2.2)*(x(1)-1);
    case 3
        y=gamma(2.8)/gamma(2.2)*(x(2)-0.5);
end
end

%benchmarking
Bench1=zeros(length(H),2);
Bench2=zeros(length(H),2);
Bench3=zeros(length(H),2);
Bench4=zeros(length(H),2);
Bench5=zeros(length(H),2);

for i=1:length(H)
    h=H(i);
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

exact = exa(t1);
Bench1(i,2)=norm((y1-exact));
Bench2(i,2)=norm((y2-exact));
Bench3(i,2)=norm((y3-exact));
Bench4(i,2)=norm((y4-exact));
Bench5(i,2)=norm((y5'-exact));
end
%%
writematrix(Bench1,'MATLAB_PIEX.csv')
writematrix(Bench2,'MATLAB_PECE.csv')
writematrix(Bench3,'MATLAB_PIRect.csv')
writematrix(Bench4,'MATLAB_PITrap.csv')
writematrix(Bench5,'MATLAB_NLFODE_VEC.csv')