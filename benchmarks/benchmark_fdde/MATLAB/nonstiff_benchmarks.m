clc

H=[2^(-3) 2^(-4) 2^(-5)  2^(-6)  2^(-7)];% 2^(-8)];% setting step-size

alpha = 0.5;
tau = 0.5;
g = @(t, y, phi) 2/gamma(3-alpha)*t^(2-alpha) - 1/gamma(2-alpha)*t^(1-alpha) + 2*tau*t - tau^2 - tau - y + phi;
u0 = 0;
phi = @(t) 0;
t0 = 0;
T = 5;

exa = @(t) t.^2-t;

%benchmarking
Bench1=zeros(length(H),2);

for i=1:length(H)
    h=H(i);
%computing the time
Bench1(i,1) = timeit(@() fdde_pi1_ex(alpha, g, tau, t0, T, phi, h));
%computing the error
[t, y] = fdde_pi1_ex(alpha, g, tau, t0, T, phi, h);

exact = exa(t);
Bench1(i,2)=norm((y-exact));
end
%%
writematrix(Bench1,'Nonstiff_MATLAB_PIEX.csv')