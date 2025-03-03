clc

H=[2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7)];% 2^(-8)];% setting step-size

alpha = [0.5, 0.5, 0.5];
beta = 0.5;
a1 = 0.5; a2 = 0.8; a3 = 1; a4 = 1; a5 = 1; a6 = 1;
sigma = [beta, 2*beta, 1+beta, 5*beta, 2, 2+beta];
A = [-10000 0 1;
     -0.05 -0.08 -0.2;
     1 0 -1];
B = [-0.6 0 0.2;
     -0.1 -0.2 0;
     0 -0.5 -0.8];
u0 = [1; 1; 1];
big_gamma = @(k) gamma(sigma(k)+1)/gamma(sigma(k)-beta+1);
g = @(t) [a1*big_gamma(1)*t^(sigma(1)-beta) + a2*big_gamma(2)*t^(sigma(2)-beta); a3*big_gamma(3)*t^(sigma(3)-beta) + a4*big_gamma(4)*t^(sigma(4)-beta); a5*big_gamma(5)*t^(sigma(5)-beta) + a6*big_gamma(6)*t^(sigma(6)-beta)] - (A + B) * [a1*t^sigma(1) + a2*t^sigma(2) + u0(1); a3*t^sigma(3) + a4*t^sigma(4) + u0(2); a5*t^sigma(5) + a6*t^sigma(6) + u0(3)];
f_fun = @(t, y) A*y + B*y + g(t);
J_fun = @(t, y) A + B;
function y=fun(t,x,k)
beta = 0.5;
u0 = [1; 1; 1];
sigma = [beta, 2*beta, 1+beta, 5*beta, 2, 2+beta];
big_gamma = @(k) gamma(sigma(k)+1)/gamma(sigma(k)-beta+1);
A = [-10000 0 1;
     -0.05 -0.08 -0.2;
     1 0 -1];
B = [-0.6 0 0.2;
     -0.1 -0.2 0;
     0 -0.5 -0.8];
a1 = 0.5; a2 = 0.8; a3 = 1; a4 = 1; a5 = 1; a6 = 1;
switch k
    case 1
        y=A(1,:)*x+B(1, :)*x + a1*big_gamma(1)*t^(sigma(1)-beta) + a2*big_gamma(2)*t^(sigma(2)-beta)- (A(1,:) + B(1,:))*[a1*t^sigma(1) + a2*t^sigma(2) + u0(1); a3*t^sigma(3) + a4*t^sigma(4) + u0(2); a5*t^sigma(5) + a6*t^sigma(6) + u0(3)];
    case 2
        y=A(2,:)*x+B(2, :)*x + a3*big_gamma(3)*t^(sigma(3)-beta) + a4*big_gamma(4)*t^(sigma(4)-beta)- (A(2,:) + B(2,:))*[a1*t^sigma(1) + a2*t^sigma(2) + u0(1); a3*t^sigma(3) + a4*t^sigma(4) + u0(2); a5*t^sigma(5) + a6*t^sigma(6) + u0(3)];
    case 3
        y=A(3,:)*x+B(3, :)*x + a5*t^sigma(5) + a6*t^sigma(6) + u0(3) - (A(3,:) + B(3,:))*[a1*t^sigma(1) + a2*t^sigma(2) + u0(1); a3*t^sigma(3) + a4*t^sigma(4) + u0(2); a5*t^sigma(5) + a6*t^sigma(6) + u0(3)];
end
end
t0 = 0;
T = 1;

exa = @(t) [a1*t.^sigma(1) + a2*t.^sigma(2) + u0(1); a3*t.^sigma(3) + a4*t.^sigma(4) + u0(2); a5*t.^sigma(5) + a6*t.^sigma(6) + u0(3)];

%benchmarking
Bench1=zeros(length(H),2);
Bench2=zeros(length(H),2);
Bench3=zeros(length(H),2);
Bench4=zeros(length(H),2);
Bench5=zeros(length(H),2);
Bench6=zeros(length(H),2);
Bench7=zeros(length(H),2);
Bench8=zeros(length(H),2);
Bench9=zeros(length(H),2);

for i=1:length(H)
    h=H(i);
%computing the time
Bench1(i,1) = timeit(@() fde_pi1_ex(alpha,f_fun,t0,T,u0,h));
Bench2(i,1) = timeit(@() fde_pi12_pc(alpha,f_fun,t0,T,u0,h));
Bench3(i,1) = timeit(@() fde_pi1_im(alpha,f_fun,J_fun,t0,T,u0,h));
Bench4(i,1) = timeit(@() fde_pi2_im(alpha,f_fun,J_fun,t0,T,u0,h));
Bench5(i,1) = timeit(@() flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],1));
Bench6(i,1) = timeit(@() flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],2));
Bench7(i,1) = timeit(@() flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],3));
Bench8(i,1) = timeit(@() nlfode_vec(@fun,alpha,u0,h,T));
Bench9(i,1) = timeit(@() pepc_nlfode(f_fun,beta,u0,h,T,1e-6,2));
%computing the error
[t1,y1]=fde_pi1_ex(alpha,f_fun,t0,T,u0,h);
[t2,y2]=fde_pi12_pc(alpha,f_fun,t0,T,u0,h);
[t3,y3]=fde_pi1_im(alpha,f_fun,J_fun,t0,T,u0,h);
[t4,y4]=fde_pi2_im(alpha,f_fun,J_fun,t0,T,u0,h);
[t5,y5]=flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],1);
[t6,y6]=flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],2);
[t7,y7]=flmm2(beta,f_fun,J_fun,t0,T,u0,h,[],3);
[y8,t8]=nlfode_vec(@fun,alpha,u0,h,T);
[y9,t9]=pepc_nlfode(f_fun,beta,u0,h,T,1e-6,2);

exact = exa(t1);
Bench1(i,2)=norm((y1-exact));
Bench2(i,2)=norm((y2-exact));
Bench3(i,2)=norm((y3-exact));
Bench4(i,2)=norm((y4-exact));
Bench5(i,2)=norm((y5-exact));
Bench6(i,2)=norm((y6-exact));
Bench7(i,2)=norm((y7-exact));
Bench8(i,2)=norm((y8'-exact));
Bench9(i,2)=norm((y9'-exact));
end
%%
writematrix(Bench1,'Stiff_MATLAB_PIEX.csv')
writematrix(Bench2,'Stiff_MATLAB_PECE.csv')
writematrix(Bench3,'Stiff_MATLAB_PIRect.csv')
writematrix(Bench4,'Stiff_MATLAB_PITrap.csv')
writematrix(Bench5,'Stiff_MATLAB_Trapzoid.csv')
writematrix(Bench6,'Stiff_MATLAB_NewtonGregory.csv')
writematrix(Bench7,'Stiff_MATLAB_BDF.csv')
% FOTF solver failed on stiff nonlinear FODE
writematrix(Bench8,'Stiff_MATLAB_NLFODE_VEC.csv')
writematrix(Bench9,'Stiff_MATLAB_PEPC_NLFODE.csv')