clc;
clear;
close all;
format long
%parameters in the model
a = 1.1;
b = 0.6;
d = 0.7;
f = 1.2;
K = 1;
D1 = 1;
D2 = 1;
m = 0.1;
c = 0.25;

L=20;
T=10;
dx = 0.2;
dt = 0.01;
x=0:dx:L;
t=0:dt:T;
Nx = length(x);
Nt = length(t);

e=dt/(dx.*dx);


%Boundary Conditions
 U(1,2:Nt) = 1;
 W(1,2:Nt) = 1;
 U(Nx-1,2:Nt) = exp(1);
 W(Nx-1,2:Nt) = exp(1);
 
 U(2,2:Nt) = U(1,2:Nt);
 W(2,2:Nt) = W(1,2:Nt);
 U(Nx,2:Nt) = U(Nx-1,2:Nt);
 W(Nx,2:Nt) = W(Nx-1,2:Nt);
 
 %Initial conditions
U(1:Nx,1) = 1.8.*(exp (-(x(1:Nx)-0.2).^2));
W(1:Nx,1) = 0.3.*(exp (-(x(1:Nx)-0.7).^2));

%---------------------------------------------------%

% Finite Difference Method

for n = 1:Nt-1
       
    for  j=3:Nx-1
    
    U(j,n+1)= e.* (U(j+1,n) - 2.*U(j,n) + U(j-1,n)) + U(j,n)+ U(j,n).*( a*dt - a*dt.*(U(j,n))./K - b*dt.*W(j,n)./(1 + m.*U(j,n))) ;
    W(j,n+1)= e.* (W(j+1,n) - 2.*W(j,n) + W(j-1,n)) + W(j,n) + W(j,n).*( - d*dt - f.*dt.*U(j,n)./(1 + m.*U(j,n)));
    
    end
end
% plot(U(:,1))
% hold on
% plot (W(:,1))
surf(U)
