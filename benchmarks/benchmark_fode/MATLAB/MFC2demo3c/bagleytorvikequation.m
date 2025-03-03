function [Y0, anan, error_norm] = bagleytorvikequation(A,B,C,D,E,G,h)
% 
% Sample function for solving Bagley-Torvik equation with zero initial conditions:
% A,B,C are coefficients of the Bagley-Torvik equation
% 
% Ay'' + By^(3/2) + Cy = F(t)
% 
% Reference: 
% I. Podlubny, Fractional Differential Equations, Academic Press, 
% San Diego, 1999, ISBN 0125588402. 
% (Correct typos: in Fig. 4 there on page 230 the results for A=B=C=1 are shown,
% replace the text "B=0.5" and "C=0.5" on page 231 with correct "B=1" and
% "C=1")

alpha1=3;
alpha2=2.5;
alpha3=2;
alpha4=1;
alpha5=0.5;

% Numerical solution:
T=0:h:100;
N=floor(100/h + 1);
M=zeros(N,N); % pre-allocate matrix for the system

% (1) Make the matrix for the alpha1-th derivative:
Dalpha1 = ban(alpha1,N,h);
Dalpha2 = ban(alpha2,N,h);
Dalpha3 = ban(alpha3,N,h);
Dalpha4 = ban(alpha4,N,h);
Dalpha5 = ban(alpha5,N,h);

% (2) Make the matrix for the second derivative:
%D2 = ban(2,N,h);

% (3) Make the matrix for the entire equation:
M=A*Dalpha1 + B.*Dalpha2 + C.*Dalpha3 + D.*Dalpha4 + E.*Dalpha5 + G.*eye(size(Dalpha1));

% Make right-hand side:
F=6*cos(T)'-(-2*T.^2-gamma(3)/(2*gamma(2.5))*T.^1.5+gamma(2)/(gamma(1.5))*T.^0.5+7)';

% (3) Utilize zero initial conditions:
M = eliminator(N,[1 2 3])*M*eliminator(N, [1 2 3])';
F= eliminator(N,[1 2 3])*F;

% (4) Solve the system BY=F:
Y=M\F;

% (5) Pre-pend the zero values (those due to zero initial conditions)
Y0 = [0; 0; 0; Y];

ana = sqrt(2)*sin(T+pi/4)+T.^2/2-T-1;
anan = ana';
error_norm = norm(anan - Y0);
disp(error_norm)
% Plot the solution:
%figure(1)
%plot(T,Y0,'g')
%grid on


