clc;
clear;
close all;
format long
%parameters in the model
b = 6.0;
a = 36.0;
g = 45.0; %The positivity of g implies that there is a constant supply of activator and inhibitor ions
sig = 9.0;
d = 1.5;

dx = 0.8;
dy = 0.8;
dt = 0.0025;
L=10;
T=1000;
x=0:dx:L;
y=0:dy:L;
t=0:dt:T;
Nx = length(x);
Nt = length(t);
Ny = length(y);

nu1=dt/(dx.*dx);
nu2=dt/(dy.*dy);

%Boundary Conditions
 U(1:Nx,1,2:Nt) = (a/5) + (g*(b-4))/(5*b);
 U(1:Nx,Ny,2:Nt) = (a/5) + (g*(b-4))/(5*b);
 U(Nx,1:Ny,2:Nt)= (a/5) + (g*(b-4))/(5*b);
 U(1,1:Ny,2:Nt)= (a/5) + (g*(b-4))/(5*b);
 V(1:Nx,1,2:Nt) = (1 + g/ b.*U(1:Nx,1,2:Nt)).*(1 + U(1:Nx,1,2:Nt).^2);
 V(1:Nx,Ny,2:Nt) = (1 + g/ b.*U(1:Nx,Ny,2:Nt)).*(1 + U(1:Nx,Ny,2:Nt).^2);
 V(Nx,1:Ny,2:Nt)= (1 + g/ b.*U(Nx,1:Ny,2:Nt)).*(1 + U(Nx,1:Ny,2:Nt).^2);
 V(1,1:Ny,2:Nt)= (1 + g/ b.*U(1,1:Ny,2:Nt)).*(1 + U(1,1:Ny,2:Nt).^2); 
 

%Initial conditions
U(1:Nx,1:Ny,1) = 0;
V(1:Nx,1:Ny,1) = 0;

% Finite Difference Method

for n = 1:Nt-1
    
    for  k=2:Nx-1
        
        for j=2:Ny-1
            U(j,k,n+1) = U(j,k,n) + a*dt - dt.*U(j,k,n) - dt.*(4.*U(j,k,n).*V(j,k,n))./(1+U(j,k,n).^2) + dt*g + nu1.*(U(j+1,k,n)-2.*U(j,k,n)+U(j-1,k,n))+ nu2.*(U(j,k+1,n)-2.*U(j,k,n)+U(j,k-1,n)); 
            V(j,k,n+1) = V(j,k,n) + sig.* ( dt*b.*(U(j,k,n)- (U(j,k,n).*V(j,k,n))./(1+U(j,k,n).^2)) + dt*g + d.*nu1.*(V(j+1,k,n)-2.*V(j,k,n)+V(j-1,k,n))+ d.*nu2.*(V(j,k+1,n)-2.*V(j,k,n)+V(j,k-1,n)));
        end
        
    end
    
end

surf(U(:,:,Nt))) %prey graph
%surf(V(:,:,Nt)) %predator graph
%contour (U) %if you want to plot contour map
